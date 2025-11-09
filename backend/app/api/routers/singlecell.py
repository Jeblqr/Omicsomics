"""Single-cell analysis API endpoints."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.pipelines.singlecell import singlecell_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(tags=["singlecell"])


# Request/Response Models
class CellRangerRequest(BaseModel):
    """Request for Cell Ranger count."""

    sample_id: int = Field(..., description="Sample ID")
    fastq_dir: str = Field(..., description="Directory containing FASTQ files")
    sample_name: str = Field(..., description="Sample identifier")
    transcriptome: str = Field(
        ..., description="Path to Cell Ranger transcriptome reference"
    )
    output_dir: str = Field(..., description="Output directory")
    chemistry: str = Field(default="auto", description="Chemistry version")
    threads: int = Field(default=16, description="Number of threads")


class ScanpyPreprocessingRequest(BaseModel):
    """Request for Scanpy preprocessing."""

    sample_id: int = Field(..., description="Sample ID")
    input_h5ad: str = Field(..., description="Input AnnData H5AD file")
    output_h5ad: str = Field(..., description="Output processed H5AD file")
    min_genes: int = Field(default=200, description="Minimum genes per cell")
    min_cells: int = Field(default=3, description="Minimum cells per gene")
    max_pct_mt: float = Field(
        default=20.0, description="Maximum mitochondrial percentage"
    )
    n_top_genes: int = Field(
        default=2000, description="Number of highly variable genes"
    )
    n_neighbors: int = Field(default=10, description="Number of neighbors for UMAP")
    leiden_resolution: float = Field(
        default=0.5, description="Leiden clustering resolution"
    )


class SeuratIntegrationRequest(BaseModel):
    """Request for Seurat batch correction."""

    sample_id: int = Field(..., description="Sample ID")
    input_h5_files: list[str] = Field(
        ..., description="List of H5 files (one per batch)"
    )
    output_rds: str = Field(..., description="Output Seurat RDS file")
    batch_key: str = Field(default="batch", description="Batch identifier column")


class CellAnnotationRequest(BaseModel):
    """Request for cell type annotation."""

    sample_id: int = Field(..., description="Sample ID")
    input_h5ad: str = Field(..., description="Input AnnData file")
    output_h5ad: str = Field(..., description="Output annotated AnnData file")
    marker_genes: dict[str, list[str]] = Field(
        ...,
        description="Dictionary mapping cell types to marker gene lists",
        example={
            "T_cells": ["CD3D", "CD3E"],
            "B_cells": ["CD19", "MS4A1"],
            "Monocytes": ["CD14", "FCGR3A"],
        },
    )


class SingleCellResponse(BaseModel):
    """Response for single-cell operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/cellranger", response_model=SingleCellResponse)
async def run_cellranger(
    request: CellRangerRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run Cell Ranger count for 10x Genomics scRNA-seq data.

    Cell Ranger performs:
    - FASTQ quality control
    - Alignment to reference transcriptome
    - UMI counting and cell barcode identification
    - Generation of feature-barcode matrices
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Cell Ranger Count",
            workflow_type="singlecell_cellranger",
            sample_id=request.sample_id,
            input_files={"fastq_dir": request.fastq_dir},
            parameters={
                "transcriptome": request.transcriptome,
                "chemistry": request.chemistry,
                "threads": request.threads,
            },
        ),
    )

    background_tasks.add_task(
        singlecell_analyzer.run_cellranger_count,
        workflow.id,
        request.fastq_dir,
        request.sample_name,
        request.transcriptome,
        request.output_dir,
        request.chemistry,
        request.threads,
        db,
    )

    return SingleCellResponse(
        workflow_id=workflow.id,
        status="queued",
        message="Cell Ranger count queued for execution",
    )


@router.post("/preprocess", response_model=SingleCellResponse)
async def preprocess_scanpy(
    request: ScanpyPreprocessingRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run Scanpy preprocessing pipeline.

    This pipeline performs:
    - Quality control (filter cells/genes)
    - Normalization and log transformation
    - Highly variable gene identification
    - Dimensionality reduction (PCA, UMAP)
    - Clustering (Leiden algorithm)
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Scanpy Preprocessing",
            workflow_type="singlecell_scanpy",
            sample_id=request.sample_id,
            input_files={"h5ad": request.input_h5ad},
            parameters={
                "min_genes": request.min_genes,
                "min_cells": request.min_cells,
                "max_pct_mt": request.max_pct_mt,
                "n_top_genes": request.n_top_genes,
                "n_neighbors": request.n_neighbors,
                "leiden_resolution": request.leiden_resolution,
            },
        ),
    )

    background_tasks.add_task(
        singlecell_analyzer.run_scanpy_preprocessing,
        workflow.id,
        request.input_h5ad,
        request.output_h5ad,
        {
            "min_genes": request.min_genes,
            "min_cells": request.min_cells,
            "max_pct_mt": request.max_pct_mt,
            "n_top_genes": request.n_top_genes,
            "n_neighbors": request.n_neighbors,
            "leiden_resolution": request.leiden_resolution,
        },
        db,
    )

    return SingleCellResponse(
        workflow_id=workflow.id,
        status="queued",
        message="Scanpy preprocessing queued for execution",
    )


@router.post("/integrate", response_model=SingleCellResponse)
async def integrate_seurat(
    request: SeuratIntegrationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run Seurat integration for batch correction.

    Seurat integration:
    - Identifies integration anchors across batches
    - Corrects batch effects while preserving biological variation
    - Enables joint analysis of multiple experiments
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Seurat Integration",
            workflow_type="singlecell_integration",
            sample_id=request.sample_id,
            input_files={"h5_files": request.input_h5_files},
            parameters={"batch_key": request.batch_key},
        ),
    )

    background_tasks.add_task(
        singlecell_analyzer.run_seurat_integration,
        workflow.id,
        request.input_h5_files,
        request.output_rds,
        request.batch_key,
        db,
    )

    return SingleCellResponse(
        workflow_id=workflow.id,
        status="queued",
        message="Seurat integration queued for execution",
    )


@router.post("/annotate", response_model=SingleCellResponse)
async def annotate_cells(
    request: CellAnnotationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Annotate cell types based on marker genes.

    This endpoint assigns cell type labels to clusters based on
    the expression of known marker genes.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Cell Type Annotation",
            workflow_type="singlecell_annotation",
            sample_id=request.sample_id,
            input_files={"h5ad": request.input_h5ad},
            parameters={"marker_genes": request.marker_genes},
        ),
    )

    background_tasks.add_task(
        singlecell_analyzer.run_cell_annotation,
        workflow.id,
        request.input_h5ad,
        request.output_h5ad,
        request.marker_genes,
        db,
    )

    return SingleCellResponse(
        workflow_id=workflow.id,
        status="queued",
        message="Cell type annotation queued for execution",
    )
