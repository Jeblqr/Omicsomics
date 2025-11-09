"""Metabolomics analysis API endpoints."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.pipelines.metabolomics import metabolomics_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(tags=["metabolomics"])


# Request/Response Models
class FeatureDetectionRequest(BaseModel):
    """Request for metabolomics feature detection."""

    sample_id: int = Field(..., description="Sample ID")
    mzml_files: list[str] = Field(..., description="mzML file paths")
    output_dir: str = Field(..., description="Output directory")
    tool: str = Field(default="xcms", description="Tool: xcms or mzmine")
    ppm: int = Field(default=25, description="Mass accuracy in ppm")
    peakwidth_min: float = Field(default=10.0, description="Minimum peak width (seconds)")
    peakwidth_max: float = Field(default=60.0, description="Maximum peak width (seconds)")


class SpectralAnnotationRequest(BaseModel):
    """Request for spectral matching and annotation."""

    sample_id: int = Field(..., description="Sample ID")
    input_file: str = Field(..., description="MGF file for GNPS or feature table for MS-DIAL")
    output_dir: str = Field(..., description="Output directory")
    tool: str = Field(default="gnps", description="Tool: gnps or msdial")
    msp_library: str | None = Field(
        default=None, description="MSP library file (required for MS-DIAL)"
    )


class QuantificationRequest(BaseModel):
    """Request for feature quantification and normalization."""

    sample_id: int = Field(..., description="Sample ID")
    feature_table: str = Field(..., description="Feature intensity table CSV")
    sample_metadata: str = Field(..., description="Sample metadata CSV")
    output_dir: str = Field(..., description="Output directory")
    normalization: str = Field(
        default="median", description="Normalization method: median, quantile, pqn"
    )


class MetabolomicsWorkflowRequest(BaseModel):
    """Complete metabolomics pipeline request."""

    sample_id: int = Field(..., description="Sample ID")
    mzml_files: list[str] = Field(..., description="mzML file paths")
    sample_metadata: str = Field(..., description="Sample metadata CSV")
    output_dir: str = Field(..., description="Output directory for all results")
    feature_tool: str = Field(default="xcms", description="Feature detection tool")
    annotation_tool: str = Field(default="gnps", description="Annotation tool")
    normalization: str = Field(default="median", description="Normalization method")
    run_annotation: bool = Field(default=True, description="Run spectral annotation")


class MetabolomicsResponse(BaseModel):
    """Response for metabolomics operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/feature-detection", response_model=MetabolomicsResponse)
async def run_feature_detection(
    request: FeatureDetectionRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run metabolomics feature detection.

    Detects metabolite features (m/z and retention time pairs) from LC-MS data.

    Supported tools:
    - XCMS: R-based tool for LC-MS data processing with advanced peak picking
    - MZmine: Java-based comprehensive platform for mass spectrometry data processing

    The process includes:
    1. Peak detection in individual samples
    2. Peak alignment across samples
    3. Retention time correction
    4. Gap filling for missing values
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Metabolomics Feature Detection ({request.tool})",
            workflow_type="metabolomics_feature_detection",
            sample_id=request.sample_id,
            input_files={"mzml_files": request.mzml_files},
            parameters={
                "tool": request.tool,
                "ppm": request.ppm,
                "peakwidth_min": request.peakwidth_min,
                "peakwidth_max": request.peakwidth_max,
            },
        ),
    )

    params = {
        "ppm": request.ppm,
        "peakwidth_min": request.peakwidth_min,
        "peakwidth_max": request.peakwidth_max,
    }

    if request.tool == "xcms":
        background_tasks.add_task(
            metabolomics_analyzer.run_xcms_feature_detection,
            workflow.id,
            request.mzml_files,
            request.output_dir,
            params,
            db,
        )
    elif request.tool == "mzmine":
        background_tasks.add_task(
            metabolomics_analyzer.run_mzmine_feature_detection,
            workflow.id,
            request.mzml_files,
            request.output_dir,
            params,
            db,
        )
    else:
        return MetabolomicsResponse(
            workflow_id=workflow.id,
            status="failed",
            message=f"Unsupported tool: {request.tool}. Use 'xcms' or 'mzmine'",
        )

    return MetabolomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.tool} feature detection queued for execution",
    )


@router.post("/spectral-annotation", response_model=MetabolomicsResponse)
async def run_spectral_annotation(
    request: SpectralAnnotationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run spectral matching and metabolite annotation.

    Annotates metabolite features by matching MS/MS spectra against spectral libraries.

    Supported tools:
    - GNPS (Global Natural Products Social Molecular Networking): 
      Web-based platform for spectral library matching with extensive natural product coverage
    - MS-DIAL: 
      Comprehensive tool for MS data processing with support for multiple spectral libraries

    The annotation helps identify known metabolites and discover novel compounds.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Metabolomics Spectral Annotation ({request.tool})",
            workflow_type="metabolomics_annotation",
            sample_id=request.sample_id,
            input_files={"input_file": request.input_file},
            parameters={
                "tool": request.tool,
                "msp_library": request.msp_library,
            },
        ),
    )

    if request.tool == "gnps":
        background_tasks.add_task(
            metabolomics_analyzer.run_gnps_annotation,
            workflow.id,
            request.input_file,
            request.output_dir,
            db,
        )
    elif request.tool == "msdial":
        if not request.msp_library:
            return MetabolomicsResponse(
                workflow_id=workflow.id,
                status="failed",
                message="MS-DIAL requires msp_library parameter",
            )
        background_tasks.add_task(
            metabolomics_analyzer.run_msdial_annotation,
            workflow.id,
            request.input_file,
            request.msp_library,
            f"{request.output_dir}/msdial_annotated.csv",
            db,
        )
    else:
        return MetabolomicsResponse(
            workflow_id=workflow.id,
            status="failed",
            message=f"Unsupported tool: {request.tool}. Use 'gnps' or 'msdial'",
        )

    return MetabolomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"{request.tool} spectral annotation queued for execution",
    )


@router.post("/quantification", response_model=MetabolomicsResponse)
async def run_quantification(
    request: QuantificationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run feature quantification and normalization.

    Processes and normalizes metabolite feature intensities across samples.

    Normalization methods:
    - median: Median normalization (adjusts each sample to have the same median intensity)
    - quantile: Quantile normalization (makes the distribution of intensities the same)
    - pqn: Probabilistic Quotient Normalization (robust to dilution differences)

    The process includes:
    1. Normalization across samples
    2. Log2 transformation
    3. Quality control metrics
    4. Export of normalized matrices
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Metabolomics Quantification ({request.normalization})",
            workflow_type="metabolomics_quantification",
            sample_id=request.sample_id,
            input_files={
                "feature_table": request.feature_table,
                "metadata": request.sample_metadata,
            },
            parameters={"normalization": request.normalization},
        ),
    )

    background_tasks.add_task(
        metabolomics_analyzer.quantify_features,
        workflow.id,
        request.feature_table,
        request.sample_metadata,
        request.output_dir,
        request.normalization,
        db,
    )

    return MetabolomicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"Feature quantification with {request.normalization} normalization queued",
    )


@router.post("/complete-pipeline", response_model=MetabolomicsResponse)
async def run_complete_metabolomics_pipeline(
    request: MetabolomicsWorkflowRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run complete metabolomics pipeline.

    This endpoint orchestrates the full metabolomics analysis workflow:
    1. Feature detection (XCMS or MZmine)
    2. Spectral annotation (GNPS or MS-DIAL) - optional
    3. Quantification and normalization

    The workflow processes untargeted LC-MS/MS data from raw mzML files
    through to normalized, annotated metabolite abundance matrices ready
    for statistical analysis.

    This comprehensive pipeline enables:
    - Discovery of differential metabolites
    - Pathway analysis
    - Biomarker identification
    - Metabolic phenotyping
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Complete Metabolomics Pipeline ({request.feature_tool})",
            workflow_type="metabolomics_complete",
            sample_id=request.sample_id,
            input_files={
                "mzml_files": request.mzml_files,
                "metadata": request.sample_metadata,
            },
            parameters={
                "feature_tool": request.feature_tool,
                "annotation_tool": request.annotation_tool,
                "normalization": request.normalization,
                "run_annotation": request.run_annotation,
            },
        ),
    )

    # Note: A complete orchestration would require sequential execution
    # For now, we just create the workflow record
    # Full implementation would use workflow engine (Nextflow/CWL)

    return MetabolomicsResponse(
        workflow_id=workflow.id,
        status="created",
        message=f"Complete metabolomics pipeline ({request.feature_tool}) created. Use workflow engine for execution.",
    )
