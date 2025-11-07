"""Genomics analysis API endpoints."""

import logging
from pathlib import Path
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.models.workflow import WorkflowStatus
from app.pipelines.genomics import genomics_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/genomics", tags=["genomics"])


# Request/Response Models
class GenomicsQCRequest(BaseModel):
    """Request for genomics QC."""

    sample_id: int = Field(..., description="Sample ID")
    fastq_files: list[str] = Field(..., description="FASTQ file paths")
    output_dir: str = Field(..., description="Output directory for QC results")


class TrimmingRequest(BaseModel):
    """Request for adapter trimming."""

    sample_id: int = Field(..., description="Sample ID")
    fastq_files: list[str] = Field(..., description="FASTQ file paths")
    output_dir: str = Field(..., description="Output directory")
    tool: str = Field(
        default="fastp", description="Trimming tool: fastp or trimmomatic"
    )
    params: dict = Field(default_factory=dict, description="Tool-specific parameters")


class AlignmentRequest(BaseModel):
    """Request for sequence alignment."""

    sample_id: int = Field(..., description="Sample ID")
    fastq_files: list[str] = Field(..., description="FASTQ file paths (trimmed)")
    reference_genome: str = Field(..., description="Path to reference genome FASTA")
    output_bam: str = Field(..., description="Output BAM file path")
    aligner: str = Field(
        default="bwa-mem", description="Aligner: bwa-mem, bowtie2, minimap2"
    )
    threads: int = Field(default=8, description="Number of threads")


class VariantCallingRequest(BaseModel):
    """Request for variant calling."""

    sample_id: int = Field(..., description="Sample ID")
    input_bam: str = Field(..., description="Input BAM file")
    reference_genome: str = Field(..., description="Reference genome FASTA")
    output_vcf: str = Field(..., description="Output VCF file path")
    caller: str = Field(
        default="gatk4", description="Variant caller: gatk4, freebayes, deepvariant"
    )


class VariantAnnotationRequest(BaseModel):
    """Request for variant annotation."""

    sample_id: int = Field(..., description="Sample ID")
    input_vcf: str = Field(..., description="Input VCF file")
    output_vcf: str = Field(..., description="Output annotated VCF")
    annotator: str = Field(default="vep", description="Annotator: vep, snpeff, annovar")


class GenomicsWorkflowRequest(BaseModel):
    """Complete genomics pipeline request."""

    sample_id: int = Field(..., description="Sample ID")
    fastq_files: list[str] = Field(..., description="Input FASTQ files")
    reference_genome: str = Field(..., description="Reference genome path")
    output_dir: str = Field(..., description="Output directory for all results")
    skip_qc: bool = Field(default=False, description="Skip FastQC step")
    skip_trimming: bool = Field(default=False, description="Skip trimming step")
    aligner: str = Field(default="bwa-mem", description="Aligner to use")
    variant_caller: str = Field(default="gatk4", description="Variant caller to use")
    annotator: str = Field(default="vep", description="Variant annotator to use")
    threads: int = Field(default=8, description="Number of threads")


class GenomicsResponse(BaseModel):
    """Response for genomics operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/qc", response_model=GenomicsResponse)
async def run_genomics_qc(
    request: GenomicsQCRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Run FastQC on genomics data."""
    # Create workflow record
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Genomics QC (FastQC)",
            workflow_type="genomics_qc",
            sample_id=request.sample_id,
            input_files={"fastq_files": request.fastq_files},
            parameters={"output_dir": request.output_dir},
        ),
    )

    # Queue background task
    background_tasks.add_task(
        genomics_analyzer.run_fastqc,
        workflow.id,
        request.fastq_files,
        request.output_dir,
        db,
    )

    return GenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="FastQC analysis queued for execution",
    )


@router.post("/trim", response_model=GenomicsResponse)
async def run_trimming(
    request: TrimmingRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Run adapter trimming on FASTQ files."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Adapter Trimming ({request.tool})",
            workflow_type="genomics_trimming",
            sample_id=request.sample_id,
            input_files={"fastq_files": request.fastq_files},
            parameters={"tool": request.tool, "params": request.params},
        ),
    )

    background_tasks.add_task(
        genomics_analyzer.run_trimming,
        workflow.id,
        request.fastq_files,
        request.output_dir,
        request.tool,
        request.params,
        db,
    )

    return GenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.tool} trimming queued for execution",
    )


@router.post("/align", response_model=GenomicsResponse)
async def run_alignment(
    request: AlignmentRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Run sequence alignment to reference genome."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Alignment ({request.aligner})",
            workflow_type="genomics_alignment",
            sample_id=request.sample_id,
            input_files={"fastq_files": request.fastq_files},
            parameters={
                "reference_genome": request.reference_genome,
                "aligner": request.aligner,
                "threads": request.threads,
            },
        ),
    )

    background_tasks.add_task(
        genomics_analyzer.run_alignment,
        workflow.id,
        request.fastq_files,
        request.reference_genome,
        request.output_bam,
        request.aligner,
        request.threads,
        db,
    )

    return GenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.aligner} alignment queued for execution",
    )


@router.post("/variant-calling", response_model=GenomicsResponse)
async def run_variant_calling(
    request: VariantCallingRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Run variant calling from aligned BAM file."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Variant Calling ({request.caller})",
            workflow_type="genomics_variant_calling",
            sample_id=request.sample_id,
            input_files={"bam": request.input_bam},
            parameters={
                "reference_genome": request.reference_genome,
                "caller": request.caller,
            },
        ),
    )

    background_tasks.add_task(
        genomics_analyzer.run_variant_calling,
        workflow.id,
        request.input_bam,
        request.reference_genome,
        request.output_vcf,
        request.caller,
        db,
    )

    return GenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.caller} variant calling queued for execution",
    )


@router.post("/annotate-variants", response_model=GenomicsResponse)
async def run_variant_annotation(
    request: VariantAnnotationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """Annotate variants with functional information."""
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Variant Annotation ({request.annotator})",
            workflow_type="genomics_annotation",
            sample_id=request.sample_id,
            input_files={"vcf": request.input_vcf},
            parameters={"annotator": request.annotator},
        ),
    )

    background_tasks.add_task(
        genomics_analyzer.run_variant_annotation,
        workflow.id,
        request.input_vcf,
        request.output_vcf,
        request.annotator,
        None,
        db,
    )

    return GenomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.annotator} annotation queued for execution",
    )


@router.post("/complete-pipeline", response_model=GenomicsResponse)
async def run_complete_genomics_pipeline(
    request: GenomicsWorkflowRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run complete genomics pipeline: QC → Trimming → Alignment → Variant Calling → Annotation.

    This endpoint orchestrates the entire genomics analysis workflow.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Complete Genomics Pipeline",
            workflow_type="genomics_complete",
            sample_id=request.sample_id,
            input_files={"fastq_files": request.fastq_files},
            parameters={
                "reference_genome": request.reference_genome,
                "skip_qc": request.skip_qc,
                "skip_trimming": request.skip_trimming,
                "aligner": request.aligner,
                "variant_caller": request.variant_caller,
                "annotator": request.annotator,
                "threads": request.threads,
            },
        ),
    )

    # Note: A complete orchestration would require sequential execution
    # For now, we just create the workflow record
    # Full implementation would use workflow engine (Nextflow/CWL)

    return GenomicsResponse(
        workflow_id=workflow.id,
        status="created",
        message="Complete genomics pipeline created. Use workflow engine for execution.",
    )
