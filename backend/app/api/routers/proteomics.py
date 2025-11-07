"""Proteomics analysis API endpoints."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.pipelines.proteomics import proteomics_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/proteomics", tags=["proteomics"])


# Request/Response Models
class RawConversionRequest(BaseModel):
    """Request for raw to mzML conversion."""

    sample_id: int = Field(..., description="Sample ID")
    raw_files: list[str] = Field(..., description="Vendor raw file paths")
    output_dir: str = Field(..., description="Output directory for mzML files")


class MaxQuantRequest(BaseModel):
    """Request for MaxQuant protein identification and quantification."""

    sample_id: int = Field(..., description="Sample ID")
    raw_files: list[str] = Field(..., description="Raw or mzML file paths")
    fasta_file: str = Field(..., description="Protein sequence database (FASTA)")
    output_dir: str = Field(..., description="Output directory")
    lfq: bool = Field(default=True, description="Enable label-free quantification")
    match_between_runs: bool = Field(default=True, description="Enable match between runs")
    max_missed_cleavages: int = Field(default=2, description="Maximum missed cleavages")


class MSFraggerRequest(BaseModel):
    """Request for MSFragger peptide identification."""

    sample_id: int = Field(..., description="Sample ID")
    mzml_files: list[str] = Field(..., description="mzML file paths")
    fasta_file: str = Field(..., description="Protein sequence database (FASTA)")
    output_dir: str = Field(..., description="Output directory")


class LFQQuantificationRequest(BaseModel):
    """Request for LFQ data processing."""

    sample_id: int = Field(..., description="Sample ID")
    protein_groups_file: str = Field(..., description="MaxQuant proteinGroups.txt file")
    output_file: str = Field(..., description="Output quantification file path")


class ProteomicsWorkflowRequest(BaseModel):
    """Complete proteomics pipeline request."""

    sample_id: int = Field(..., description="Sample ID")
    raw_files: list[str] = Field(..., description="Vendor raw file paths")
    fasta_file: str = Field(..., description="Protein sequence database (FASTA)")
    output_dir: str = Field(..., description="Output directory for all results")
    tool: str = Field(default="maxquant", description="Tool: maxquant or msfragger")
    lfq: bool = Field(default=True, description="Enable label-free quantification")
    convert_raw: bool = Field(default=True, description="Convert raw files to mzML first")


class ProteomicsResponse(BaseModel):
    """Response for proteomics operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/convert-raw", response_model=ProteomicsResponse)
async def convert_raw_to_mzml(
    request: RawConversionRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Convert vendor raw files to mzML format.

    Supports Thermo raw files using ThermoRawFileParser.
    The mzML format is an open standard for mass spectrometry data.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Raw to mzML Conversion",
            workflow_type="proteomics_conversion",
            sample_id=request.sample_id,
            input_files={"raw_files": request.raw_files},
            parameters={"output_dir": request.output_dir},
        ),
    )

    background_tasks.add_task(
        proteomics_analyzer.convert_raw_to_mzml,
        workflow.id,
        request.raw_files,
        request.output_dir,
        db,
    )

    return ProteomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="Raw to mzML conversion queued for execution",
    )


@router.post("/maxquant", response_model=ProteomicsResponse)
async def run_maxquant_analysis(
    request: MaxQuantRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run MaxQuant for protein identification and quantification.

    MaxQuant is a widely-used tool for analyzing mass spectrometry data.
    It performs peptide/protein identification and label-free quantification (LFQ).

    Features:
    - Peptide and protein identification
    - Label-free quantification (LFQ)
    - Match between runs for increased coverage
    - Post-translational modification analysis
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="MaxQuant Analysis",
            workflow_type="proteomics_maxquant",
            sample_id=request.sample_id,
            input_files={"raw_files": request.raw_files, "fasta": request.fasta_file},
            parameters={
                "lfq": request.lfq,
                "match_between_runs": request.match_between_runs,
                "max_missed_cleavages": request.max_missed_cleavages,
            },
        ),
    )

    background_tasks.add_task(
        proteomics_analyzer.run_maxquant,
        workflow.id,
        request.raw_files,
        request.fasta_file,
        request.output_dir,
        {
            "lfq": request.lfq,
            "match_between_runs": request.match_between_runs,
            "max_missed_cleavages": request.max_missed_cleavages,
        },
        db,
    )

    return ProteomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="MaxQuant analysis queued for execution",
    )


@router.post("/msfragger", response_model=ProteomicsResponse)
async def run_msfragger_analysis(
    request: MSFraggerRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run MSFragger for peptide identification.

    MSFragger is a fast database search tool for proteomics.
    It's particularly useful for open searches and identifying post-translational modifications.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="MSFragger Analysis",
            workflow_type="proteomics_msfragger",
            sample_id=request.sample_id,
            input_files={"mzml_files": request.mzml_files, "fasta": request.fasta_file},
            parameters={},
        ),
    )

    background_tasks.add_task(
        proteomics_analyzer.run_msfragGer,
        workflow.id,
        request.mzml_files,
        request.fasta_file,
        request.output_dir,
        db,
    )

    return ProteomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="MSFragger analysis queued for execution",
    )


@router.post("/lfq-quantification", response_model=ProteomicsResponse)
async def process_lfq_quantification(
    request: LFQQuantificationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Process label-free quantification (LFQ) data.

    Takes MaxQuant proteinGroups.txt output and processes LFQ intensities:
    - Filters contaminants and reverse hits
    - Extracts LFQ intensity columns
    - Applies log2 transformation
    - Generates clean quantification matrix
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="LFQ Quantification Processing",
            workflow_type="proteomics_lfq",
            sample_id=request.sample_id,
            input_files={"protein_groups": request.protein_groups_file},
            parameters={},
        ),
    )

    background_tasks.add_task(
        proteomics_analyzer.quantify_lfq,
        workflow.id,
        request.protein_groups_file,
        request.output_file,
        db,
    )

    return ProteomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message="LFQ quantification processing queued",
    )


@router.post("/complete-pipeline", response_model=ProteomicsResponse)
async def run_complete_proteomics_pipeline(
    request: ProteomicsWorkflowRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run complete proteomics pipeline.

    This endpoint orchestrates the full proteomics analysis workflow:
    1. Raw file conversion (optional)
    2. Peptide/protein identification (MaxQuant or MSFragger)
    3. Quantification processing (if LFQ enabled)

    The workflow handles mass spectrometry data from raw files through
    to quantified protein abundances.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Complete Proteomics Pipeline ({request.tool})",
            workflow_type="proteomics_complete",
            sample_id=request.sample_id,
            input_files={"raw_files": request.raw_files, "fasta": request.fasta_file},
            parameters={
                "tool": request.tool,
                "lfq": request.lfq,
                "convert_raw": request.convert_raw,
            },
        ),
    )

    # Note: A complete orchestration would require sequential execution
    # For now, we just create the workflow record
    # Full implementation would use workflow engine (Nextflow/CWL)

    return ProteomicsResponse(
        workflow_id=workflow.id,
        status="created",
        message=f"Complete proteomics pipeline ({request.tool}) created. Use workflow engine for execution.",
    )
