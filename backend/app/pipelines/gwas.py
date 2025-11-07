"""Genome-Wide Association Study (GWAS) analysis pipeline."""

import asyncio
import logging
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class GWASAnalyzer:
    """GWAS analysis pipeline orchestrator."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize GWAS analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_gwas")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def run_plink_qc(
        self,
        workflow_id: int,
        bed_file: str,
        bim_file: str,
        fam_file: str,
        output_prefix: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run PLINK quality control on genotype data.

        Args:
            workflow_id: Workflow database ID
            bed_file: PLINK BED file
            bim_file: PLINK BIM file
            fam_file: PLINK FAM file
            output_prefix: Output file prefix
            params: QC parameters
            db: Database session

        Returns:
            Result with QC-filtered files
        """
        params = params or {}
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

        # QC parameters
        geno = params.get("geno", 0.02)  # SNP missing rate
        mind = params.get("mind", 0.02)  # Individual missing rate
        maf = params.get("maf", 0.01)  # Minor allele frequency
        hwe = params.get("hwe", 1e-6)  # Hardy-Weinberg equilibrium

        cmd = [
            "plink",
            "--bfile",
            bed_file.replace(".bed", ""),
            "--geno",
            str(geno),
            "--mind",
            str(mind),
            "--maf",
            str(maf),
            "--hwe",
            str(hwe),
            "--make-bed",
            "--out",
            output_prefix,
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running PLINK QC: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "bed": f"{output_prefix}.bed",
                    "bim": f"{output_prefix}.bim",
                    "fam": f"{output_prefix}.fam",
                    "log": f"{output_prefix}.log",
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
                            error_message=f"PLINK QC failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"PLINK QC failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_association_test(
        self,
        workflow_id: int,
        bed_file: str,
        phenotype_file: str,
        output_prefix: str,
        covariates_file: str | None = None,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run GWAS association test using PLINK.

        Args:
            workflow_id: Workflow database ID
            bed_file: PLINK BED file (QC-filtered)
            phenotype_file: Phenotype file
            output_prefix: Output file prefix
            covariates_file: Optional covariates file
            params: Association test parameters
            db: Database session

        Returns:
            Result with association results
        """
        params = params or {}
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "plink",
            "--bfile",
            bed_file.replace(".bed", ""),
            "--pheno",
            phenotype_file,
            "--assoc",
            "--adjust",
            "--out",
            output_prefix,
        ]

        # Add covariates if provided
        if covariates_file:
            cmd.extend(["--covar", covariates_file])

        # Linear vs logistic regression
        if params.get("binary_trait", False):
            cmd.append("--logistic")
        else:
            cmd.append("--linear")

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running GWAS: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "assoc": f"{output_prefix}.assoc"
                    + (".logistic" if params.get("binary_trait") else ".linear"),
                    "adjusted": f"{output_prefix}.assoc.adjusted",
                    "log": f"{output_prefix}.log",
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
                            error_message=f"GWAS failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"GWAS failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def calculate_ld(
        self,
        workflow_id: int,
        bed_file: str,
        output_prefix: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Calculate linkage disequilibrium (LD) matrix.

        Args:
            workflow_id: Workflow database ID
            bed_file: PLINK BED file
            output_prefix: Output file prefix
            params: LD parameters
            db: Database session

        Returns:
            Result with LD matrix
        """
        params = params or {}
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

        # LD parameters
        ld_window = params.get("ld_window", 1000)  # Window size in kb
        ld_window_r2 = params.get("ld_window_r2", 0.2)  # r2 threshold

        cmd = [
            "plink",
            "--bfile",
            bed_file.replace(".bed", ""),
            "--r2",
            "--ld-window",
            str(ld_window),
            "--ld-window-r2",
            str(ld_window_r2),
            "--out",
            output_prefix,
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Calculating LD: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "ld": f"{output_prefix}.ld",
                    "log": f"{output_prefix}.log",
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
                            error_message=f"LD calculation failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"LD calculation failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_prs_calculation(
        self,
        workflow_id: int,
        bed_file: str,
        weights_file: str,
        output_prefix: str,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Calculate Polygenic Risk Score (PRS).

        Args:
            workflow_id: Workflow database ID
            bed_file: PLINK BED file (target genotypes)
            weights_file: SNP weights file (from GWAS summary stats)
            output_prefix: Output file prefix
            db: Database session

        Returns:
            Result with PRS scores
        """
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "plink",
            "--bfile",
            bed_file.replace(".bed", ""),
            "--score",
            weights_file,
            "--out",
            output_prefix,
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Calculating PRS: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "profile": f"{output_prefix}.profile",
                    "log": f"{output_prefix}.log",
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
                            error_message=f"PRS calculation failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"PRS calculation failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_mtag_analysis(
        self,
        workflow_id: int,
        summary_stats_files: dict[str, str],
        output_dir: str,
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run Multi-Trait Analysis of GWAS (MTAG) for cross-trait analysis.

        Args:
            workflow_id: Workflow database ID
            summary_stats_files: Dictionary mapping trait names to summary stats files
            output_dir: Output directory
            params: MTAG parameters
            db: Database session

        Returns:
            Result with MTAG results
        """
        params = params or {}
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Generate MTAG config file
        config_content = f"""
# MTAG Configuration
output_dir: {output_dir}
traits:
"""
        for trait_name, sumstats_file in summary_stats_files.items():
            config_content += f"  - name: {trait_name}\n"
            config_content += f"    file: {sumstats_file}\n"

        config_file = self.work_dir / f"workflow_{workflow_id}_mtag_config.yaml"
        config_file.write_text(config_content)

        # Python script for MTAG
        python_script = f"""
import pandas as pd
import numpy as np
from scipy import stats
import json

# Load summary statistics
traits = {list(summary_stats_files.keys())}
sumstats = {{}}

for trait in traits:
    sumstats[trait] = pd.read_csv("{summary_stats_files[trait]}", sep='\\t')
    print(f"Loaded {{trait}}: {{len(sumstats[trait])}} SNPs")

# Simple cross-trait meta-analysis (placeholder for MTAG)
# Full MTAG implementation requires complex genetic correlation estimation

# Find common SNPs across all traits
common_snps = set(sumstats[traits[0]]['SNP'])
for trait in traits[1:]:
    common_snps &= set(sumstats[trait]['SNP'])

print(f"Common SNPs across all traits: {{len(common_snps)}}")

# Perform inverse-variance weighted meta-analysis
results = []
for snp in common_snps:
    snp_data = {{}}
    for trait in traits:
        trait_data = sumstats[trait][sumstats[trait]['SNP'] == snp].iloc[0]
        snp_data[trait] = {{
            'beta': trait_data.get('BETA', trait_data.get('B', 0)),
            'se': trait_data.get('SE', 1),
            'p': trait_data.get('P', 1)
        }}
    
    # Inverse variance weighted meta-analysis
    weights = [1 / (snp_data[t]['se'] ** 2) for t in traits]
    meta_beta = sum(snp_data[t]['beta'] * w for t, w in zip(traits, weights)) / sum(weights)
    meta_se = 1 / np.sqrt(sum(weights))
    meta_z = meta_beta / meta_se
    meta_p = 2 * (1 - stats.norm.cdf(abs(meta_z)))
    
    results.append({{
        'SNP': snp,
        'meta_beta': meta_beta,
        'meta_se': meta_se,
        'meta_z': meta_z,
        'meta_p': meta_p
    }})

# Save results
results_df = pd.DataFrame(results)
results_df = results_df.sort_values('meta_p')
results_df.to_csv("{output_dir}/mtag_results.csv", index=False)

# Summary
summary = {{
    'n_traits': len(traits),
    'n_common_snps': len(common_snps),
    'n_significant': sum(results_df['meta_p'] < 5e-8),
    'top_snp': results_df.iloc[0].to_dict() if len(results_df) > 0 else None
}}

with open("{output_dir}/mtag_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("MTAG analysis completed")
print(f"Significant SNPs (p < 5e-8): {{summary['n_significant']}}")
"""

        script_path = self.work_dir / f"workflow_{workflow_id}_mtag.py"
        script_path.write_text(python_script)

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            cmd = ["python", str(script_path)]
            logger.info(f"Running MTAG: {' '.join(cmd)}")

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "results": f"{output_dir}/mtag_results.csv",
                    "summary": f"{output_dir}/mtag_summary.json",
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
                            error_message=f"MTAG failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"MTAG failed: {e}")
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
gwas_analyzer = GWASAnalyzer()
