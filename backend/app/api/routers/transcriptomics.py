"""Transcriptomics analysis API endpoints."""

import logging
from pathlib import Path
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.models.workflow import WorkflowStatus
from app.pipelines.transcriptomics import transcriptomics_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(tags=["transcriptomics"])


# Request/Response Models
class QuantificationRequest(BaseModel):
    """Request for RNA-seq quantification."""

    sample_id: int = Field(..., description="Sample ID")
    fastq_files: list[str] = Field(..., description="FASTQ file paths")
    reference_index: str = Field(..., description="Path to reference index")
    output_dir: str = Field(..., description="Output directory")
    tool: str = Field(
        default="salmon", description="Tool: salmon, kallisto, star, hisat2"
    )
    threads: int = Field(default=8, description="Number of threads")


class CountMatrixRequest(BaseModel):
    """Request for count matrix generation."""

    sample_id: int = Field(..., description="Sample ID")
    input_bam: str = Field(..., description="Input BAM file")
    gtf_file: str = Field(..., description="Gene annotation GTF file")
    output_file: str = Field(..., description="Output counts file")
    threads: int = Field(default=4, description="Number of threads")


class DifferentialExpressionRequest(BaseModel):
    """Request for differential expression analysis."""

    sample_id: int = Field(..., description="Sample ID")
    counts_matrix: str = Field(..., description="Count matrix file path")
    sample_metadata: str = Field(..., description="Sample metadata CSV path")
    output_dir: str = Field(..., description="Output directory")
    tool: str = Field(default="deseq2", description="Tool: deseq2, edger, limma")
    design_formula: str = Field(default="~ condition", description="R design formula")


class TranscriptomicsResponse(BaseModel):
    """Response for transcriptomics operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/quantify", response_model=TranscriptomicsResponse)
async def run_quantification(
    request: QuantificationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run RNA-seq quantification.

    Supports multiple tools:
    - salmon: Fast, alignment-free quantification
    - kallisto: Fast pseudo-alignment
    - star: Gold-standard aligner with quantification
    - hisat2: Alignment-based approach
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"RNA-seq Quantification ({request.tool})",
            workflow_type="transcriptomics_quantification",
            sample_id=request.sample_id,
            input_files={"fastq_files": request.fastq_files},
            parameters={
                "reference_index": request.reference_index,
                "tool": request.tool,
                "threads": request.threads,
            },
        ),
    )

    background_tasks.add_task(
        transcriptomics_analyzer.run_alignment_quantification,
        workflow.id,
        request.fastq_files,
        request.reference_index,
        request.output_dir,
        request.tool,
        request.threads,
        db,
    )

    return TranscriptomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.tool} quantification queued for execution",
    )


@router.post("/count-matrix", response_model=TranscriptomicsResponse)
async def generate_count_matrix(
    request: CountMatrixRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Generate count matrix from aligned BAM file using featureCounts.

    This endpoint takes an aligned BAM file and a GTF annotation file
    to generate a gene-level count matrix suitable for differential expression analysis.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Count Matrix Generation (featureCounts)",
            workflow_type="transcriptomics_counts",
            sample_id=request.sample_id,
            input_files={"bam": request.input_bam, "gtf": request.gtf_file},
            parameters={"threads": request.threads},
        ),
    )

    background_tasks.add_task(
        transcriptomics_analyzer.run_feature_counts,
        workflow.id,
        request.input_bam,
        request.gtf_file,
        request.output_file,
        request.threads,
        db,
    )

    return TranscriptomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="featureCounts queued for execution",
    )


@router.post("/differential-expression", response_model=TranscriptomicsResponse)
async def run_differential_expression(
    request: DifferentialExpressionRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run differential expression analysis.

    This endpoint performs statistical analysis to identify differentially expressed genes
    between conditions. It uses DESeq2 by default and generates:
    - Results table with log2 fold changes and adjusted p-values
    - MA plot
    - PCA plot
    - Volcano plot

    The sample_metadata CSV should have a 'condition' column (or match your design formula).
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Differential Expression ({request.tool})",
            workflow_type="transcriptomics_de",
            sample_id=request.sample_id,
            input_files={
                "counts_matrix": request.counts_matrix,
                "metadata": request.sample_metadata,
            },
            parameters={"tool": request.tool, "design_formula": request.design_formula},
        ),
    )

    background_tasks.add_task(
        transcriptomics_analyzer.run_differential_expression,
        workflow.id,
        request.counts_matrix,
        request.sample_metadata,
        request.output_dir,
        request.tool,
        request.design_formula,
        db,
    )

    return TranscriptomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.tool} differential expression analysis queued",
    )


@router.get("/enrichment/gsea")
async def run_gsea_enrichment(
    gene_list: str,
    organism: str = "human",
    database: str = "GO_Biological_Process_2021",
    current_user: Annotated[User, Depends(get_current_user)] = None,
):
    """
    Run Gene Set Enrichment Analysis (GSEA).

    This endpoint is a placeholder for GSEA integration.
    Full implementation would use GSEApy or call Enrichr API.
    """
    return {
        "status": "not_implemented",
        "message": "GSEA enrichment analysis not yet implemented",
        "suggestion": "Use GSEApy Python package or Enrichr web API",
    }
