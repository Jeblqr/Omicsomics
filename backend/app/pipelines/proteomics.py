"""Proteomics analysis pipeline for mass spectrometry data."""

import asyncio
import logging
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class ProteomicsAnalyzer:
    """Proteomics analysis pipeline orchestrator for MS data."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize proteomics analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_proteomics")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def convert_raw_to_mzml(
        self,
        workflow_id: int,
        raw_files: list[str],
        output_dir: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Convert vendor raw files to mzML format.

        Args:
            workflow_id: Workflow database ID
            raw_files: List of vendor raw files
            output_dir: Output directory for mzML files
            db: Database session

        Returns:
            Result with mzML file paths
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            mzml_files = []
            logs = ""

            for raw_file in raw_files:
                raw_path = Path(raw_file)
                output_mzml = Path(output_dir) / f"{raw_path.stem}.mzML"
                
                # ThermoRawFileParser for Thermo files
                cmd = [
                    "ThermoRawFileParser",
                    "-i", raw_file,
                    "-o", output_dir,
                    "-f", "2",  # mzML format
                ]

                logger.info(f"Converting {raw_file} to mzML: {' '.join(cmd)}")
                process = await asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )

                stdout, stderr = await process.communicate()
                logs += f"\n=== {raw_file} ===\n"
                logs += stdout.decode() + "\n" + stderr.decode()

                if process.returncode == 0:
                    mzml_files.append(str(output_mzml))
                else:
                    if db:
                        await workflow_service.update_workflow(
                            db,
                            workflow_id,
                            workflow_schema.WorkflowUpdate(
                                status=WorkflowStatus.FAILED,
                                logs=logs,
                                error_message=f"Conversion failed for {raw_file}",
                            ),
                        )
                    return {"status": "failed", "error": f"Conversion failed for {raw_file}"}

            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.COMPLETED,
                        logs=logs,
                        output_files={"mzml_files": mzml_files},
                    ),
                )
            return {"status": "success", "mzml_files": mzml_files}

        except Exception as e:
            logger.error(f"Raw to mzML conversion failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_maxquant(
        self,
        workflow_id: int,
        raw_files: list[str],
        fasta_file: str,
        output_dir: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run MaxQuant for peptide/protein identification and quantification.

        Args:
            workflow_id: Workflow database ID
            raw_files: List of raw/mzML files
            fasta_file: Protein sequence database (FASTA)
            output_dir: Output directory
            params: MaxQuant parameters
            db: Database session

        Returns:
            Result with MaxQuant output files
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        params = params or {}

        # Generate mqpar.xml configuration
        mqpar_xml = self._generate_maxquant_config(
            raw_files, fasta_file, output_dir, params
        )
        
        mqpar_path = Path(output_dir) / "mqpar.xml"
        mqpar_path.write_text(mqpar_xml)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            # Run MaxQuant
            cmd = ["maxquant", str(mqpar_path)]

            logger.info(f"Running MaxQuant: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=output_dir,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "protein_groups": f"{output_dir}/combined/txt/proteinGroups.txt",
                    "peptides": f"{output_dir}/combined/txt/peptides.txt",
                    "evidence": f"{output_dir}/combined/txt/evidence.txt",
                    "summary": f"{output_dir}/combined/txt/summary.txt",
                }

                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files=output_files,
                        ),
                    )
                return {"status": "success", **output_files}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"MaxQuant failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"MaxQuant failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    def _generate_maxquant_config(
        self,
        raw_files: list[str],
        fasta_file: str,
        output_dir: str,
        params: dict[str, Any],
    ) -> str:
        """Generate MaxQuant mqpar.xml configuration file."""
        # Simplified mqpar.xml template
        file_list = "\n".join([f"<string>{f}</string>" for f in raw_files])
        
        xml = f"""<?xml version="1.0" encoding="utf-8"?>
<MaxQuantParams xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <filePaths>
    {file_list}
  </filePaths>
  <fastaFiles>
    <string>{fasta_file}</string>
  </fastaFiles>
  <fixedModifications>
    <string>Carbamidomethyl (C)</string>
  </fixedModifications>
  <variableModifications>
    <string>Oxidation (M)</string>
    <string>Acetyl (Protein N-term)</string>
  </variableModifications>
  <enzymes>
    <string>Trypsin/P</string>
  </enzymes>
  <maxMissedCleavages>{params.get('max_missed_cleavages', 2)}</maxMissedCleavages>
  <lfqMode>{1 if params.get('lfq', True) else 0}</lfqMode>
  <matchBetweenRuns>{1 if params.get('match_between_runs', True) else 0}</matchBetweenRuns>
</MaxQuantParams>"""
        return xml

    async def run_msfragGer(
        self,
        workflow_id: int,
        mzml_files: list[str],
        fasta_file: str,
        output_dir: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run MSFragger for peptide identification.

        Args:
            workflow_id: Workflow database ID
            mzml_files: List of mzML files
            fasta_file: Protein sequence database
            output_dir: Output directory
            db: Database session

        Returns:
            Result with MSFragger output
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            # MSFragger command
            cmd = [
                "java", "-jar", "MSFragger.jar",
                fasta_file,
                *mzml_files,
            ]

            logger.info(f"Running MSFragger: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=output_dir,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                pepxml_files = [f"{Path(f).stem}.pepXML" for f in mzml_files]
                
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={"pepxml_files": pepxml_files},
                        ),
                    )
                return {"status": "success", "pepxml_files": pepxml_files}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"MSFragger failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"MSFragger failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def quantify_lfq(
        self,
        workflow_id: int,
        protein_groups_file: str,
        output_file: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Process LFQ (Label-Free Quantification) data.

        Args:
            workflow_id: Workflow database ID
            protein_groups_file: MaxQuant proteinGroups.txt
            output_file: Output quantification file
            db: Database session

        Returns:
            Result with quantification data
        """
        # Python script for LFQ processing
        python_script = f"""
import pandas as pd
import numpy as np

# Load protein groups
df = pd.read_csv("{protein_groups_file}", sep='\\t')

# Filter contaminants and reverse hits
df = df[df['Potential contaminant'].isna()]
df = df[df['Reverse'].isna()]

# Extract LFQ intensities
lfq_cols = [col for col in df.columns if col.startswith('LFQ intensity')]
lfq_data = df[['Protein IDs', 'Gene names'] + lfq_cols].copy()

# Log2 transform
for col in lfq_cols:
    lfq_data[col] = np.log2(lfq_data[col].replace(0, np.nan))

# Save
lfq_data.to_csv("{output_file}", index=False)

print(f"Processed {{len(lfq_data)}} protein groups")
print(f"LFQ intensity columns: {{len(lfq_cols)}}")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_lfq.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Processing LFQ: {' '.join(cmd)}")
            
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={"lfq_quantification": output_file},
                        ),
                    )
                return {"status": "success", "lfq_quantification": output_file}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"LFQ processing failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"LFQ processing failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}


# Global instance
proteomics_analyzer = ProteomicsAnalyzer()
