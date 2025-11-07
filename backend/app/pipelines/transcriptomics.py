"""Transcriptomics analysis pipeline for bulk RNA-seq."""

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


class TranscriptomicsAnalyzer:
    """Transcriptomics analysis pipeline orchestrator for bulk RNA-seq."""

    def __init__(self, work_dir: Path | None = None):
        """Initialize transcriptomics analyzer."""
        self.work_dir = work_dir or Path("/tmp/omicsomics_transcriptomics")
        self.work_dir.mkdir(parents=True, exist_ok=True)

    async def run_alignment_quantification(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_index: str,
        output_dir: str,
        tool: str = "salmon",
        threads: int = 8,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run alignment and quantification for RNA-seq.

        Args:
            workflow_id: Workflow database ID
            input_files: FASTQ files (R1, R2 for paired-end)
            reference_index: Path to reference index (STAR/Salmon/Kallisto)
            output_dir: Output directory
            tool: 'star', 'hisat2', 'salmon', or 'kallisto'
            threads: Number of threads
            db: Database session

        Returns:
            Result with quantification files
        """
        if tool == "salmon":
            return await self._run_salmon(
                workflow_id, input_files, reference_index, output_dir, threads, db
            )
        elif tool == "kallisto":
            return await self._run_kallisto(
                workflow_id, input_files, reference_index, output_dir, threads, db
            )
        elif tool == "star":
            return await self._run_star(
                workflow_id, input_files, reference_index, output_dir, threads, db
            )
        elif tool == "hisat2":
            return await self._run_hisat2(
                workflow_id, input_files, reference_index, output_dir, threads, db
            )
        else:
            return {"status": "error", "error": f"Unknown tool: {tool}"}

    async def _run_salmon(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_index: str,
        output_dir: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Salmon for transcript quantification."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if len(input_files) == 1:
            # Single-end
            cmd = [
                "salmon",
                "quant",
                "-i",
                reference_index,
                "-l",
                "A",  # Auto-detect library type
                "-r",
                input_files[0],
                "-o",
                output_dir,
                "-p",
                str(threads),
            ]
        elif len(input_files) == 2:
            # Paired-end
            cmd = [
                "salmon",
                "quant",
                "-i",
                reference_index,
                "-l",
                "A",
                "-1",
                input_files[0],
                "-2",
                input_files[1],
                "-o",
                output_dir,
                "-p",
                str(threads),
            ]
        else:
            return {"status": "error", "error": "Invalid number of input files"}

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running Salmon: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                quant_file = str(Path(output_dir) / "quant.sf")
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={
                                "quant_file": quant_file,
                                "output_dir": output_dir,
                            },
                        ),
                    )
                return {
                    "status": "success",
                    "quant_file": quant_file,
                    "output_dir": output_dir,
                }
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Salmon failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Salmon failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_kallisto(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_index: str,
        output_dir: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run Kallisto for transcript quantification."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if len(input_files) == 1:
            # Single-end requires fragment length and SD
            cmd = [
                "kallisto",
                "quant",
                "-i",
                reference_index,
                "-o",
                output_dir,
                "-t",
                str(threads),
                "--single",
                "-l",
                "200",  # Average fragment length
                "-s",
                "20",  # Fragment length SD
                input_files[0],
            ]
        elif len(input_files) == 2:
            # Paired-end
            cmd = [
                "kallisto",
                "quant",
                "-i",
                reference_index,
                "-o",
                output_dir,
                "-t",
                str(threads),
                input_files[0],
                input_files[1],
            ]
        else:
            return {"status": "error", "error": "Invalid number of input files"}

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running Kallisto: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                abundance_file = str(Path(output_dir) / "abundance.tsv")
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={
                                "abundance_file": abundance_file,
                                "output_dir": output_dir,
                            },
                        ),
                    )
                return {"status": "success", "abundance_file": abundance_file}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"Kallisto failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"Kallisto failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_star(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_index: str,
        output_dir: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run STAR aligner for RNA-seq."""
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        cmd = [
            "STAR",
            "--runThreadN",
            str(threads),
            "--genomeDir",
            reference_index,
            "--readFilesIn",
            *input_files,
            "--outFileNamePrefix",
            f"{output_dir}/",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--quantMode",
            "GeneCounts",
        ]

        # Handle gzipped files
        if any(f.endswith(".gz") for f in input_files):
            cmd.extend(["--readFilesCommand", "zcat"])

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running STAR: {' '.join(cmd)}")
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            stdout, stderr = await process.communicate()
            logs = stdout.decode() + "\n" + stderr.decode()

            if process.returncode == 0:
                bam_file = f"{output_dir}/Aligned.sortedByCoord.out.bam"
                counts_file = f"{output_dir}/ReadsPerGene.out.tab"
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.COMPLETED,
                            logs=logs,
                            output_files={
                                "bam": bam_file,
                                "counts": counts_file,
                                "output_dir": output_dir,
                            },
                        ),
                    )
                return {"status": "success", "bam": bam_file, "counts": counts_file}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"STAR failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"STAR failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def _run_hisat2(
        self,
        workflow_id: int,
        input_files: list[str],
        reference_index: str,
        output_dir: str,
        threads: int,
        db: AsyncSession,
    ) -> dict[str, Any]:
        """Run HISAT2 for RNA-seq alignment."""
        # Placeholder for HISAT2 implementation
        return {"status": "error", "error": "HISAT2 not yet implemented"}

    async def run_feature_counts(
        self,
        workflow_id: int,
        input_bam: str,
        gtf_file: str,
        output_file: str,
        threads: int = 4,
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Generate count matrix using featureCounts.

        Args:
            workflow_id: Workflow database ID
            input_bam: Input BAM file(s)
            gtf_file: Gene annotation GTF file
            output_file: Output counts file
            threads: Number of threads
            db: Database session

        Returns:
            Result with counts file path
        """
        cmd = [
            "featureCounts",
            "-T",
            str(threads),
            "-a",
            gtf_file,
            "-o",
            output_file,
            input_bam,
        ]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running featureCounts: {' '.join(cmd)}")
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
                            output_files={"counts": output_file},
                        ),
                    )
                return {"status": "success", "counts": output_file}
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"featureCounts failed: exit code {process.returncode}",
                        ),
                    )
                return {"status": "failed", "error": f"Exit code {process.returncode}"}

        except Exception as e:
            logger.error(f"featureCounts failed: {e}")
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(
                        status=WorkflowStatus.FAILED, error_message=str(e)
                    ),
                )
            return {"status": "error", "error": str(e)}

    async def run_differential_expression(
        self,
        workflow_id: int,
        counts_matrix: str,
        sample_metadata: str,
        output_dir: str,
        tool: str = "deseq2",
        design_formula: str = "~ condition",
        db: AsyncSession | None = None,
    ) -> dict[str, Any]:
        """
        Run differential expression analysis.

        Args:
            workflow_id: Workflow database ID
            counts_matrix: Count matrix file
            sample_metadata: Sample metadata CSV
            output_dir: Output directory
            tool: 'deseq2', 'edger', or 'limma'
            design_formula: R design formula
            db: Database session

        Returns:
            Result with DE results files
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # R script for DESeq2 analysis
        r_script = f"""
library(DESeq2)
library(ggplot2)

# Load data
counts <- read.table("{counts_matrix}", header=TRUE, row.names=1)
coldata <- read.csv("{sample_metadata}", row.names=1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = {design_formula})

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), file="{output_dir}/deseq2_results.csv")

# MA plot
pdf("{output_dir}/MA_plot.pdf")
plotMA(res)
dev.off()

# PCA plot
vsd <- vst(dds, blind=FALSE)
pdf("{output_dir}/PCA_plot.pdf")
plotPCA(vsd, intgroup="condition")
dev.off()

# Volcano plot
pdf("{output_dir}/volcano_plot.pdf")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano Plot"))
with(subset(res, padj<.05), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
"""

        script_path = Path(output_dir) / "run_deseq2.R"
        script_path.write_text(r_script)

        cmd = ["Rscript", str(script_path)]

        try:
            if db:
                await workflow_service.update_workflow(
                    db,
                    workflow_id,
                    workflow_schema.WorkflowUpdate(status=WorkflowStatus.RUNNING),
                )

            logger.info(f"Running DESeq2: {' '.join(cmd)}")
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
                                "results": f"{output_dir}/deseq2_results.csv",
                                "ma_plot": f"{output_dir}/MA_plot.pdf",
                                "pca_plot": f"{output_dir}/PCA_plot.pdf",
                                "volcano_plot": f"{output_dir}/volcano_plot.pdf",
                            },
                        ),
                    )
                return {
                    "status": "success",
                    "results": f"{output_dir}/deseq2_results.csv",
                    "plots": {
                        "ma": f"{output_dir}/MA_plot.pdf",
                        "pca": f"{output_dir}/PCA_plot.pdf",
                        "volcano": f"{output_dir}/volcano_plot.pdf",
                    },
                }
            else:
                if db:
                    await workflow_service.update_workflow(
                        db,
                        workflow_id,
                        workflow_schema.WorkflowUpdate(
                            status=WorkflowStatus.FAILED,
                            logs=logs,
                            error_message=f"DESeq2 failed: exit code {process.returncode}",
                        ),
                    )
                return {
                    "status": "failed",
                    "error": f"Exit code {process.returncode}",
                    "logs": logs,
                }

        except Exception as e:
            logger.error(f"DESeq2 failed: {e}")
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
transcriptomics_analyzer = TranscriptomicsAnalyzer()
