"""Epigenomics analysis API endpoints for ChIP-seq, ATAC-seq, and methylation data."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.pipelines.epigenomics import epigenomics_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(tags=["epigenomics"])


# Request/Response Models
class EpigenomicsAlignmentRequest(BaseModel):
    """Request for ChIP-seq/ATAC-seq alignment."""

    sample_id: int = Field(..., description="Sample ID")
    input_files: list[str] = Field(..., description="FASTQ file paths (R1, R2 for paired-end)")
    reference_genome: str = Field(..., description="Path to reference genome index")
    output_bam: str = Field(..., description="Output BAM file path")
    aligner: str = Field(default="bowtie2", description="Aligner: bowtie2 or bwa")
    threads: int = Field(default=8, description="Number of threads")


class PeakCallingRequest(BaseModel):
    """Request for peak calling."""

    sample_id: int = Field(..., description="Sample ID")
    treatment_bam: str = Field(..., description="Treatment BAM file")
    control_bam: str | None = Field(None, description="Control/Input BAM file (optional)")
    output_dir: str = Field(..., description="Output directory")
    peak_caller: str = Field(default="macs2", description="Peak caller: macs2 or macs3")
    genome_size: str = Field(default="hs", description="Genome size: 'hs' (human), 'mm' (mouse), or integer")
    broad: bool = Field(default=False, description="Call broad peaks")
    q_value: float | None = Field(None, description="Q-value cutoff")
    p_value: float | None = Field(None, description="P-value cutoff")
    nomodel: bool = Field(default=False, description="Bypass model building")
    shift: int | None = Field(None, description="Shift size")
    extsize: int | None = Field(None, description="Extension size")


class MotifAnalysisRequest(BaseModel):
    """Request for motif analysis."""

    sample_id: int = Field(..., description="Sample ID")
    peak_file: str = Field(..., description="Peak file (BED or narrowPeak format)")
    genome_fasta: str = Field(..., description="Reference genome FASTA file")
    output_dir: str = Field(..., description="Output directory")
    tool: str = Field(default="homer", description="Tool: homer or meme")
    motif_length: list[int] = Field(default=[8, 10, 12], description="Motif lengths to search")


class BigWigRequest(BaseModel):
    """Request for BigWig generation."""

    sample_id: int = Field(..., description="Sample ID")
    input_bam: str = Field(..., description="Input BAM file")
    output_bigwig: str = Field(..., description="Output BigWig file path")
    genome_sizes: str = Field(..., description="Genome chromosome sizes file")
    normalize: bool = Field(default=True, description="Normalize using RPKM")


class EpigenomicsWorkflowRequest(BaseModel):
    """Complete epigenomics pipeline request."""

    sample_id: int = Field(..., description="Sample ID")
    treatment_fastq: list[str] = Field(..., description="Treatment FASTQ files")
    control_fastq: list[str] | None = Field(None, description="Control/Input FASTQ files")
    reference_genome: str = Field(..., description="Reference genome index path")
    genome_fasta: str = Field(..., description="Reference genome FASTA")
    genome_sizes: str = Field(..., description="Genome chromosome sizes file")
    output_dir: str = Field(..., description="Output directory for all results")
    aligner: str = Field(default="bowtie2", description="Aligner to use")
    peak_caller: str = Field(default="macs2", description="Peak caller to use")
    genome_size: str = Field(default="hs", description="Genome size for MACS2")
    run_motif_analysis: bool = Field(default=False, description="Run motif analysis")
    generate_bigwig: bool = Field(default=True, description="Generate BigWig files")
    threads: int = Field(default=8, description="Number of threads")


class EpigenomicsResponse(BaseModel):
    """Response for epigenomics operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/align", response_model=EpigenomicsResponse)
async def run_epigenomics_alignment(
    request: EpigenomicsAlignmentRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Align ChIP-seq/ATAC-seq reads to reference genome."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Epigenomics Alignment ({request.aligner})",
            workflow_type="epigenomics_alignment",
            sample_id=request.sample_id,
            input_files={"fastq_files": request.input_files},
            parameters={
                "reference_genome": request.reference_genome,
                "aligner": request.aligner,
                "threads": request.threads,
            },
        ),
    )

    background_tasks.add_task(
        epigenomics_analyzer.run_alignment,
        workflow.id,
        request.input_files,
        request.reference_genome,
        request.output_bam,
        request.aligner,
        request.threads,
        db,
    )

    return EpigenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.aligner} alignment queued for execution",
    )


@router.post("/peak-calling", response_model=EpigenomicsResponse)
async def run_peak_calling(
    request: PeakCallingRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Call peaks from ChIP-seq/ATAC-seq data using MACS2/MACS3."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Peak Calling ({request.peak_caller})",
            workflow_type="epigenomics_peak_calling",
            sample_id=request.sample_id,
            input_files={
                "treatment_bam": request.treatment_bam,
                "control_bam": request.control_bam,
            },
            parameters={
                "peak_caller": request.peak_caller,
                "genome_size": request.genome_size,
                "broad": request.broad,
                "q_value": request.q_value,
                "p_value": request.p_value,
                "nomodel": request.nomodel,
                "shift": request.shift,
                "extsize": request.extsize,
            },
        ),
    )

    # Prepare parameters
    params = {}
    if request.broad:
        params["broad"] = True
    if request.q_value is not None:
        params["q_value"] = request.q_value
    if request.p_value is not None:
        params["p_value"] = request.p_value
    if request.nomodel:
        params["nomodel"] = True
    if request.shift is not None:
        params["shift"] = request.shift
    if request.extsize is not None:
        params["extsize"] = request.extsize

    background_tasks.add_task(
        epigenomics_analyzer.run_peak_calling,
        workflow.id,
        request.treatment_bam,
        request.control_bam,
        request.output_dir,
        request.peak_caller,
        request.genome_size,
        params,
        db,
    )

    return EpigenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.peak_caller} peak calling queued for execution",
    )


@router.post("/motif-analysis", response_model=EpigenomicsResponse)
async def run_motif_analysis(
    request: MotifAnalysisRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Perform motif analysis on peak regions using HOMER or MEME."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Motif Analysis ({request.tool})",
            workflow_type="epigenomics_motif",
            sample_id=request.sample_id,
            input_files={"peak_file": request.peak_file},
            parameters={
                "tool": request.tool,
                "motif_length": request.motif_length,
            },
        ),
    )

    background_tasks.add_task(
        epigenomics_analyzer.run_motif_analysis,
        workflow.id,
        request.peak_file,
        request.genome_fasta,
        request.output_dir,
        request.tool,
        request.motif_length,
        db,
    )

    return EpigenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.tool} motif analysis queued for execution",
    )


@router.post("/bigwig", response_model=EpigenomicsResponse)
async def generate_bigwig(
    request: BigWigRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Generate BigWig file for genome browser visualization."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="BigWig Generation",
            workflow_type="epigenomics_bigwig",
            sample_id=request.sample_id,
            input_files={"bam": request.input_bam},
            parameters={
                "normalize": request.normalize,
            },
        ),
    )

    background_tasks.add_task(
        epigenomics_analyzer.generate_bigwig,
        workflow.id,
        request.input_bam,
        request.output_bigwig,
        request.genome_sizes,
        request.normalize,
        db,
    )

    return EpigenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="BigWig generation queued for execution",
    )


@router.post("/complete-pipeline", response_model=EpigenomicsResponse)
async def run_complete_epigenomics_pipeline(
    request: EpigenomicsWorkflowRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run complete epigenomics pipeline: Alignment → Peak Calling → Motif Analysis → BigWig Generation.

    This endpoint orchestrates the entire epigenomics analysis workflow for ChIP-seq/ATAC-seq data.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Complete Epigenomics Pipeline",
            workflow_type="epigenomics_complete",
            sample_id=request.sample_id,
            input_files={
                "treatment_fastq": request.treatment_fastq,
                "control_fastq": request.control_fastq,
            },
            parameters={
                "reference_genome": request.reference_genome,
                "genome_fasta": request.genome_fasta,
                "genome_sizes": request.genome_sizes,
                "aligner": request.aligner,
                "peak_caller": request.peak_caller,
                "genome_size": request.genome_size,
                "run_motif_analysis": request.run_motif_analysis,
                "generate_bigwig": request.generate_bigwig,
                "threads": request.threads,
            },
        ),
    )

    # Note: A complete orchestration would require sequential execution
    # Full implementation would use workflow engine (Nextflow/CWL) to coordinate steps
    # For now, we create the workflow record that can be executed by the workflow engine

    return EpigenomicsResponse(
        workflow_id=workflow.id,
        status="created",
        message="Complete epigenomics pipeline created. Use workflow engine for execution.",
    )
