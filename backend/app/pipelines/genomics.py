"""Genomics analysis pipeline for WGS/WES data."""

import asyncio
import json
import logging
from pathlib import Path
from typing import Any

from app.models.workflow import WorkflowStatus
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)


class GenomicsAnalyzer:
    """Genomics analysis pipeline orchestrator."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize genomics analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_genomics")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def run_fastqc(
        self,
        workflow_id: int,
        input_files: list[str],
        output_dir: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """
        Run FastQC quality control on FASTQ files.

        Args:
            workflow_id: Workflow database ID
            input_files: List of FASTQ file paths
            output_dir: Output directory
            db: Database session

        Returns:
            Result dictionary with status and output paths
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        cmd = ["fastqc", "-o", output_dir, "-t", "4"] + input_files

        try:
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
            )

            logger.info(f"Running FastQC: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.COMPLETED,
                        logs=logs,
                        output_files={"qc_dir": output_dir},
                    ),
                )
                return {"status": "success", "qc_dir": output_dir}
            else:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED,
                        logs=logs,
                        error_message=f"FastQC failed: exit code {process.returncode}",
                    ),
                )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"FastQC failed: {e}")
            await workflow_service.update_workflow(
                db,
                workflow_id,
                workflow_schema.WorkflowUpdate(
                    status=WorkflowStatus.FAILED, error_message=str(e)
                ),
            )
            return {"status": "error", "error": str(e)}

    async def run_trimming(
        self,
        workflow_id: int,
        input_files: list[str],
        output_dir: str,
        tool: str = "fastp",
        params: dict[str, Any] | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run adapter trimming with fastp or Trimmomatic.

        Args:
            workflow_id: Workflow database ID
            input_files: FASTQ files (R1, R2 for paired-end)
            output_dir: Output directory
            tool: 'fastp' or 'trimmomatic'
            params: Tool-specific parameters
            db: Database session

        Returns:
            Result with trimmed file paths
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        params = params or {}

        if tool == "fastp":
            return await self._run_fastp(
                workflow_id, input_files, output_dir, params, db
            )
        elif tool == "trimmomatic":
            return await self._run_trimmomatic(
                workflow_id, input_files, output_dir, params, db
            )
        else:
            return {"status": "error", "error": f"Unknown trimming tool: {tool}"}

    async def _run_fastp(
        self,
        workflow_id: int,
        input_files: list[str],
        output_dir: str,
        params: dict[str, Any],
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run fastp for adapter trimming and QC."""
        output_prefix = Path(output_dir) / "trimmed"

        if len(input_files) == 1:
            # Single-end
            cmd = [
                "fastp",
                "-i",
                input_files[0],
                "-o",
                f"{output_prefix}_R1.fastq.gz",
                "-j",
                f"{output_prefix}_fastp.json",
                "-h",
                f"{output_prefix}_fastp.html",
            ]
        elif len(input_files) == 2:
            # Paired-end
            cmd = [
                "fastp",
                "-i",
                input_files[0],
                "-I",
                input_files[1],
                "-o",
                f"{output_prefix}_R1.fastq.gz",
                "-O",
                f"{output_prefix}_R2.fastq.gz",
                "-j",
                f"{output_prefix}_fastp.json",
                "-h",
                f"{output_prefix}_fastp.html",
            ]
        else:
            return {"status": "error", "error": "Invalid number of input files"}

        # Add custom parameters
        for key, value in params.items():
            cmd.extend([f"--{key}", str(value)])

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running fastp: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                output_files = {
                    "trimmed_r1": f"{output_prefix}_R1.fastq.gz",
                    "report_json": f"{output_prefix}_fastp.json",
                    "report_html": f"{output_prefix}_fastp.html",
                }
                if len(input_files) == 2:
                    output_files["trimmed_r2"] = f"{output_prefix}_R2.fastq.gz"

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
                            error_message=f"fastp failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"fastp failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_trimmomatic(
        self,
        workflow_id: int,
        input_files: list[str],
        output_dir: str,
        params: dict[str, Any],
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Trimmomatic for adapter trimming."""
        # Placeholder for Trimmomatic implementation
        return {"status": "error", "error": "Trimmomatic not yet implemented"}

    async def run_alignment(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        aligner: str = "bwa-mem",
        threads: int = 8,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Align reads to reference genome.

        Args:
            workflow_id: Workflow database ID
            input_files: FASTQ files (R1, R2 for paired-end)
            reference_genome: Path to reference genome FASTA (indexed)
            output_bam: Output BAM file path
            aligner: 'bwa-mem', 'bowtie2', or 'minimap2'
            threads: Number of threads
            db: Database session

        Returns:
            Result with BAM file path
        """
        if aligner == "bwa-mem":
            return await self._run_bwa_mem(
                workflow_id, input_files, reference_genome, output_bam, threads, db
            )
        elif aligner == "bowtie2":
            return await self._run_bowtie2(
                workflow_id, input_files, reference_genome, output_bam, threads, db
            )
        elif aligner == "minimap2":
            return await self._run_minimap2(
                workflow_id, input_files, reference_genome, output_bam, threads, db
            )
        else:
            return {"status": "error", "error": f"Unknown aligner: {aligner}"}

    async def _run_bwa_mem(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run BWA-MEM alignment."""
        Path(output_bam).parent.mkdir(parents=True, exist_ok=True)

        # BWA MEM alignment with piping to samtools for BAM conversion
        if len(input_files) == 1:
            cmd = f"bwa mem -t {threads} {reference_genome} {input_files[0]} | samtools sort -@ {threads} -o {output_bam}"
        elif len(input_files) == 2:
            cmd = f"bwa mem -t {threads} {reference_genome} {input_files[0]} {input_files[1]} | samtools sort -@ {threads} -o {output_bam}"
        else:
            return {"status": "error", "error": "Invalid number of input files"}

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running BWA-MEM: {cmd}")
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
                            output_files={
                                "bam": output_bam,
                                "bai": f"{output_bam}.bai",
                            },
                        ),
                    )
                return {
                    "status": "success",
                    "bam": output_bam,
                    "bai": f"{output_bam}.bai",
                }
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"BWA-MEM failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"BWA-MEM failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_bowtie2(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Bowtie2 alignment."""
        # Placeholder for Bowtie2 implementation
        return {"status": "error", "error": "Bowtie2 not yet implemented"}

    async def _run_minimap2(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_genome: str,
        output_bam: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Minimap2 alignment."""
        # Placeholder for Minimap2 implementation
        return {"status": "error", "error": "Minimap2 not yet implemented"}

    async def run_variant_calling(
        self,
        workflow_id: int,
        input_bam: str,
        reference_genome: str,
        output_vcf: str,
        caller: str = "gatk4",
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Call variants from aligned BAM file.

        Args:
            workflow_id: Workflow database ID
            input_bam: Input BAM file
            reference_genome: Reference genome FASTA
            output_vcf: Output VCF file path
            caller: 'gatk4', 'freebayes', or 'deepvariant'
            db: Database session

        Returns:
            Result with VCF file path
        """
        if caller == "gatk4":
            return await self._run_gatk_haplotypecaller(
                workflow_id, input_bam, reference_genome, output_vcf, db
            )
        elif caller == "freebayes":
            return await self._run_freebayes(
                workflow_id, input_bam, reference_genome, output_vcf, db
            )
        elif caller == "deepvariant":
            return await self._run_deepvariant(
                workflow_id, input_bam, reference_genome, output_vcf, db
            )
        else:
            return {"status": "error", "error": f"Unknown variant caller: {caller}"}

    async def _run_gatk_haplotypecaller(
        self,
        workflow_id: int,
        input_bam: str,
        reference_genome: str,
        output_vcf: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run GATK4 HaplotypeCaller."""
        Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            "gatk",
            "HaplotypeCaller",
            "-R",
            reference_genome,
            "-I",
            input_bam,
            "-O",
            output_vcf,
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running GATK HaplotypeCaller: {' '.join(cmd)}")
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
                            output_files={"vcf": output_vcf},
                        ),
                    )
                return {"status": "success", "vcf": output_vcf}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"GATK failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"GATK HaplotypeCaller failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_freebayes(
        self,
        workflow_id: int,
        input_bam: str,
        reference_genome: str,
        output_vcf: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run FreeBayes variant calling."""
        # Placeholder for FreeBayes implementation
        return {"status": "error", "error": "FreeBayes not yet implemented"}

    async def _run_deepvariant(
        self,
        workflow_id: int,
        input_bam: str,
        reference_genome: str,
        output_vcf: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run DeepVariant variant calling."""
        # Placeholder for DeepVariant implementation
        return {"status": "error", "error": "DeepVariant not yet implemented"}

    async def run_variant_annotation(
        self,
        workflow_id: int,
        input_vcf: str,
        output_vcf: str,
        annotator: str = "vep",
        reference_genome: str | None = None,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Annotate variants with functional information.

        Args:
            workflow_id: Workflow database ID
            input_vcf: Input VCF file
            output_vcf: Output annotated VCF
            annotator: 'vep', 'snpeff', or 'annovar'
            reference_genome: Reference genome (if needed)
            db: Database session

        Returns:
            Result with annotated VCF path
        """
        if annotator == "vep":
            return await self._run_vep(workflow_id, input_vcf, output_vcf, db)
        elif annotator == "snpeff":
            return await self._run_snpeff(workflow_id, input_vcf, output_vcf, db)
        elif annotator == "annovar":
            return await self._run_annovar(workflow_id, input_vcf, output_vcf, db)
        else:
            return {"status": "error", "error": f"Unknown annotator: {annotator}"}

    async def _run_vep(
        self,
        workflow_id: int,
        input_vcf: str,
        output_vcf: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Ensembl VEP annotation."""
        cmd = [
            "vep",
            "-i",
            input_vcf,
            "-o",
            output_vcf,
            "--vcf",
            "--cache",
            "--everything",
            "--fork",
            "4",
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running VEP: {' '.join(cmd)}")
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
                            output_files={"annotated_vcf": output_vcf},
                        ),
                    )
                return {"status": "success", "annotated_vcf": output_vcf}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"VEP failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"VEP failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_snpeff(
        self,
        workflow_id: int,
        input_vcf: str,
        output_vcf: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run SnpEff annotation."""
        # Placeholder for SnpEff implementation
        return {"status": "error", "error": "SnpEff not yet implemented"}

    async def _run_annovar(
        self,
        workflow_id: int,
        input_vcf: str,
        output_vcf: str,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run ANNOVAR annotation."""
        # Placeholder for ANNOVAR implementation
        return {"status": "error", "error": "ANNOVAR not yet implemented"}


# Global instance
genomics_analyzer = GenomicsAnalyzer()
