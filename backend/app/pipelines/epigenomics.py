"""Epigenomics analysis pipeline for ChIP-seq, ATAC-seq, and methylation data."""

import asyncio
import logging
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class EpigenomicsAnalyzer:
    """Epigenomics analysis pipeline orchestrator."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize epigenomics analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_epigenomics")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def run_alignment(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        aligner: str = "bowtie2",
        threads: int = 8,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Align ChIP-seq/ATAC-seq reads to reference genome.

        Args:
            workflow_id: Workflow database ID
            input_files: FASTQ files (R1, R2 for paired-end)
            reference_genome: Path to reference genome index
            output_bam: Output BAM file path
            aligner: 'bowtie2' or 'bwa'
            threads: Number of threads
            db: Database session

        Returns:
            Result with BAM file path
        """
        if aligner == "bowtie2":
            return await self._run_bowtie2(
                workflow_id, input_files, reference_genome, output_bam, threads, db
            )
        elif aligner == "bwa":
            return await self._run_bwa(
                workflow_id, input_files, reference_genome, output_bam, threads, db
            )
        else:
            return {"status": "error", "error": f"Unknown aligner: {aligner}"}

    async def _run_bowtie2(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Bowtie2 alignment for ChIP-seq/ATAC-seq."""
        Path(output_bam).parent.mkdir(parents=True, exist_ok=True)

        if len(input_files) == 1:
            # Single-end
            cmd = f"bowtie2 -p {threads} -x {reference_genome} -U {input_files[0]} | samtools sort -@ {threads} -o {output_bam}"
        elif len(input_files) == 2:
            # Paired-end
            cmd = f"bowtie2 -p {threads} -x {reference_genome} -1 {input_files[0]} -2 {input_files[1]} | samtools sort -@ {threads} -o {output_bam}"
        else:
            return {"status": "error", "error": "Invalid number of input files"}

        try:
            if db:
                await workflow_service.update_workflow(
                    db, workflow_id, workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING)
                )

            logger.info(f"Running Bowtie2: {cmd}")
            process = await asyncio.create_subprocess_shell(
                cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                # Index BAM file
                index_cmd = ["samtools", "index", output_bam]
                await asyncio.create_subprocess_exec(*index_cmd)

                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={"bam": output_bam, "bai": f"{output_bam}.bai"},
                        ),
                    )
                return {"status": "success", "bam": output_bam, "bai": f"{output_bam}.bai"}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Bowtie2 failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Bowtie2 failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_bwa(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run BWA alignment for ChIP-seq/ATAC-seq."""
        # Similar to BWA-MEM implementation in genomics
        return {"status": "error", "error": "BWA for epigenomics not yet implemented"}

    async def run_peak_calling(
        self,
        workflow_id: int,
        treatment_bam: str,
        control_bam: str | None,
        output_dir: str,
        peak_caller: str = "macs2",
        genome_size: str = "hs",
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Call peaks from ChIP-seq/ATAC-seq data.

        Args:
            workflow_id: Workflow database ID
            treatment_bam: Treatment BAM file
            control_bam: Control/Input BAM file (optional)
            output_dir: Output directory
            peak_caller: 'macs2' or 'macs3'
            genome_size: Genome size ('hs' for human, 'mm' for mouse, or integer)
            params: Additional parameters
            db: Database session

        Returns:
            Result with peak files
        """
        if peak_caller in ["macs2", "macs3"]:
            return await self._run_macs(
                workflow_id, treatment_bam, control_bam, output_dir, genome_size, params or {}, db
            )
        else:
            return {"status": "error", "error": f"Unknown peak caller: {peak_caller}"}

    async def _run_macs(
        self,
        workflow_id: int,
        treatment_bam: str,
        control_bam: str | None,
        output_dir: str,
        genome_size: str,
        params: dict[str, Any],
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run MACS2/3 for peak calling."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Build MACS2 command
        cmd = [
            "macs2", "callpeak",
            "-t", treatment_bam,
            "-f", "BAM",
            "-g", genome_size,
            "-n", "peaks",
            "--outdir", output_dir,
        ]

        # Add control if provided
        if control_bam:
            cmd.extend(["-c", control_bam])

        # Add additional parameters
        if params.get("broad"):
            cmd.append("--broad")
        if params.get("q_value"):
            cmd.extend(["-q", str(params["q_value"])])
        if params.get("p_value"):
            cmd.extend(["-p", str(params["p_value"])])
        if params.get("nomodel"):
            cmd.append("--nomodel")
        if params.get("shift"):
            cmd.extend(["--shift", str(params["shift"])])
        if params.get("extsize"):
            cmd.extend(["--extsize", str(params["extsize"])])

        try:
            if db:
                await workflow_service.update_workflow(
                    db, workflow_id, workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING)
                )

            logger.info(f"Running MACS2: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "peaks": f"{output_dir}/peaks_peaks.narrowPeak",
                    "summits": f"{output_dir}/peaks_summits.bed",
                    "xls": f"{output_dir}/peaks_peaks.xls",
                }
                if params.get("broad"):
                    output_files["broad_peaks"] = f"{output_dir}/peaks_peaks.broadPeak"

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
                            error_message=f"MACS2 failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}", "logs": logs}

        except Exception as e:
            logger.error(f"MACS2 failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_motif_analysis(
        self,
        workflow_id: int,
        peak_file: str,
        genome_fasta: str,
        output_dir: str,
        tool: str = "homer",
        motif_length: list[int] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Perform motif analysis on peak regions.

        Args:
            workflow_id: Workflow database ID
            peak_file: Peak file (BED or narrowPeak format)
            genome_fasta: Reference genome FASTA
            output_dir: Output directory
            tool: 'homer' or 'meme'
            motif_length: List of motif lengths to search
            db: Database session

        Returns:
            Result with motif analysis results
        """
        if tool == "homer":
            return await self._run_homer(
                workflow_id, peak_file, genome_fasta, output_dir, motif_length or [8, 10, 12], db
            )
        elif tool == "meme":
            return await self._run_meme(
                workflow_id, peak_file, genome_fasta, output_dir, db
            )
        else:
            return {"status": "error", "error": f"Unknown motif analysis tool: {tool}"}

    async def _run_homer(
        self,
        workflow_id: int,
        peak_file: str,
        genome_fasta: str,
        output_dir: str,
        motif_length: list[int],
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run HOMER for motif analysis."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # HOMER findMotifsGenome.pl command
        motif_len_str = ",".join(map(str, motif_length))
        cmd = [
            "findMotifsGenome.pl",
            peak_file,
            genome_fasta,
            output_dir,
            "-size", "200",
            "-len", motif_len_str,
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db, workflow_id, workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING)
                )

            logger.info(f"Running HOMER: {' '.join(cmd)}")
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
                            output_files={
                                "motif_html": f"{output_dir}/homerResults.html",
                                "known_motifs": f"{output_dir}/knownResults.html",
                                "output_dir": output_dir,
                            },
                        ),
                    )
                return {"status": "success", "output_dir": output_dir}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"HOMER failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}", "logs": logs}

        except Exception as e:
            logger.error(f"HOMER failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_meme(
        self,
        workflow_id: int,
        peak_file: str,
        genome_fasta: str,
        output_dir: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run MEME Suite for motif analysis."""
        # Placeholder for MEME implementation
        return {"status": "error", "error": "MEME not yet implemented"}

    async def generate_bigwig(
        self,
        workflow_id: int,
        input_bam: str,
        output_bigwig: str,
        genome_sizes: str,
        normalize: bool = True,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Generate BigWig file for visualization.

        Args:
            workflow_id: Workflow database ID
            input_bam: Input BAM file
            output_bigwig: Output BigWig file path
            genome_sizes: Genome chromosome sizes file
            normalize: Whether to normalize (RPKM)
            db: Database session

        Returns:
            Result with BigWig file path
        """
        try:
            if db:
                await workflow_service.update_workflow(
                    db, workflow_id, workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING)
                )

            # bamCoverage from deepTools
            cmd = [
                "bamCoverage",
                "-b", input_bam,
                "-o", output_bigwig,
                "--binSize", "10",
            ]

            if normalize:
                cmd.extend(["--normalizeUsing", "RPKM"])

            logger.info(f"Generating BigWig: {' '.join(cmd)}")
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
                            output_files={"bigwig": output_bigwig},
                        ),
                    )
                return {"status": "success", "bigwig": output_bigwig}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"bamCoverage failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"BigWig generation failed: {e}")
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
epigenomics_analyzer = EpigenomicsAnalyzer()
