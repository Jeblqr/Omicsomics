"""Multi-omics integration API endpoints."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.pipelines.multiomics import multiomics_integrator
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/multiomics", tags=["multiomics"])


# Request/Response Models
class MOFA2Request(BaseModel):
    """Request for MOFA2 multi-omics integration."""

    sample_id: int = Field(..., description="Sample ID")
    data_matrices: dict[str, str] = Field(
        ...,
        description="Dictionary mapping omics type to data matrix CSV file path",
        example={
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
            "metabolomics": "/data/metabolite_features.csv",
        },
    )
    output_dir: str = Field(..., description="Output directory")
    n_factors: int = Field(default=10, description="Number of latent factors")
    convergence_mode: str = Field(
        default="fast", description="Convergence mode: fast, medium, or slow"
    )


class DIABLORequest(BaseModel):
    """Request for DIABLO multi-omics integration."""

    sample_id: int = Field(..., description="Sample ID")
    data_matrices: dict[str, str] = Field(
        ...,
        description="Dictionary mapping omics type to data matrix CSV file path",
    )
    phenotype_file: str = Field(
        ..., description="CSV file with sample phenotype/class labels"
    )
    output_dir: str = Field(..., description="Output directory")
    n_components: int = Field(default=2, description="Number of latent components")
    design_correlation: float = Field(
        default=0.1, description="Design matrix correlation value (0-1)"
    )


class PathwayEnrichmentRequest(BaseModel):
    """Request for multi-omics pathway enrichment."""

    sample_id: int = Field(..., description="Sample ID")
    feature_lists: dict[str, str] = Field(
        ...,
        description="Dictionary mapping omics type to feature list CSV file",
        example={
            "transcriptomics": "/data/deg_list.csv",
            "proteomics": "/data/dap_list.csv",
        },
    )
    organism: str = Field(
        default="hsapiens", description="Organism: hsapiens, mmusculus, etc."
    )
    output_dir: str = Field(..., description="Output directory")


class SampleMatchingRequest(BaseModel):
    """Request for matching samples across omics datasets."""

    sample_id: int = Field(..., description="Sample ID")
    omics_tables: dict[str, str] = Field(
        ..., description="Dictionary mapping omics type to data table CSV"
    )
    output_file: str = Field(..., description="Output sample mapping JSON file")


class MultiOmicsWorkflowRequest(BaseModel):
    """Complete multi-omics integration workflow request."""

    sample_id: int = Field(..., description="Sample ID")
    data_matrices: dict[str, str] = Field(..., description="Omics data matrices")
    phenotype_file: str | None = Field(
        default=None, description="Phenotype file (required for DIABLO)"
    )
    output_dir: str = Field(..., description="Output directory for all results")
    integration_method: str = Field(
        default="mofa2", description="Integration method: mofa2 or diablo"
    )
    run_pathway_enrichment: bool = Field(
        default=True, description="Run pathway enrichment analysis"
    )
    organism: str = Field(
        default="hsapiens", description="Organism for pathway analysis"
    )


class MultiOmicsResponse(BaseModel):
    """Response for multi-omics operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/mofa2", response_model=MultiOmicsResponse)
async def run_mofa2_integration(
    request: MOFA2Request,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run MOFA2 (Multi-Omics Factor Analysis v2) for unsupervised data integration.

    MOFA2 identifies latent factors that capture shared and unique variation
    across multiple omics layers. It's particularly useful for:
    - Discovering molecular patterns across omics
    - Identifying drivers of heterogeneity
    - Dimensionality reduction of multi-omics data
    - Imputing missing values across datasets

    The method decomposes each omics matrix into a set of latent factors,
    where each factor has:
    - Factor values: sample scores on the factor
    - Weights: feature loadings indicating importance
    - Variance explained: contribution to each omics layer

    MOFA2 handles:
    - Different sample sizes across omics
    - Missing values
    - Different data distributions
    - Multiple sample groups
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="MOFA2 Multi-Omics Integration",
            workflow_type="multiomics_mofa2",
            sample_id=request.sample_id,
            input_files={"data_matrices": request.data_matrices},
            parameters={
                "n_factors": request.n_factors,
                "convergence_mode": request.convergence_mode,
            },
        ),
    )

    params = {
        "n_factors": request.n_factors,
        "convergence_mode": request.convergence_mode,
    }

    background_tasks.add_task(
        multiomics_integrator.run_mofa2_integration,
        workflow.id,
        request.data_matrices,
        request.output_dir,
        params,
        db,
    )

    return MultiOmicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"MOFA2 integration with {len(request.data_matrices)} omics layers queued",
    )


@router.post("/diablo", response_model=MultiOmicsResponse)
async def run_diablo_integration(
    request: DIABLORequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run DIABLO for supervised multi-omics integration and biomarker discovery.

    DIABLO (Data Integration Analysis for Biomarker discovery using Latent cOmponents)
    is a supervised method that integrates multiple omics while maximizing
    discrimination between sample groups (e.g., disease vs. control).

    Key features:
    - Supervised: uses phenotype/class labels
    - Variable selection: identifies key biomarkers in each omics
    - Maximizes covariance between omics and phenotype
    - Handles high-dimensional data with feature selection

    The method produces:
    - Latent components that discriminate classes
    - Selected features (biomarkers) for each omics
    - Correlation networks across omics
    - Classification performance metrics

    Use cases:
    - Disease biomarker discovery
    - Stratification of patient groups
    - Understanding multi-omics signatures
    - Predictive modeling with multi-omics
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="DIABLO Multi-Omics Integration",
            workflow_type="multiomics_diablo",
            sample_id=request.sample_id,
            input_files={
                "data_matrices": request.data_matrices,
                "phenotype": request.phenotype_file,
            },
            parameters={
                "n_components": request.n_components,
                "design_correlation": request.design_correlation,
            },
        ),
    )

    params = {
        "n_components": request.n_components,
        "design_correlation": request.design_correlation,
    }

    background_tasks.add_task(
        multiomics_integrator.run_diablo_integration,
        workflow.id,
        request.data_matrices,
        request.phenotype_file,
        request.output_dir,
        params,
        db,
    )

    return MultiOmicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"DIABLO integration with {len(request.data_matrices)} omics layers queued",
    )


@router.post("/pathway-enrichment", response_model=MultiOmicsResponse)
async def run_pathway_enrichment(
    request: PathwayEnrichmentRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run pathway enrichment analysis across multiple omics layers.

    Performs enrichment analysis on feature lists from different omics types
    to identify enriched biological pathways, GO terms, and other functional
    categories.

    This multi-omics approach reveals:
    - Shared pathways across omics
    - Omics-specific biological processes
    - Integrated biological insights
    - System-level understanding

    Supported enrichment sources:
    - KEGG pathways
    - GO (Gene Ontology) terms
    - Reactome pathways
    - WikiPathways
    - Custom gene sets

    The analysis helps interpret multi-omics results by connecting
    molecular features to biological function.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Multi-Omics Pathway Enrichment",
            workflow_type="multiomics_pathways",
            sample_id=request.sample_id,
            input_files={"feature_lists": request.feature_lists},
            parameters={"organism": request.organism},
        ),
    )

    background_tasks.add_task(
        multiomics_integrator.run_pathway_enrichment,
        workflow.id,
        request.feature_lists,
        request.organism,
        request.output_dir,
        db,
    )

    return MultiOmicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"Pathway enrichment for {len(request.feature_lists)} omics layers queued",
    )


@router.post("/match-samples", response_model=MultiOmicsResponse)
async def match_samples(
    request: SampleMatchingRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Match samples across different omics datasets.

    Before integrating multi-omics data, it's essential to identify which
    samples are measured across all omics layers. This endpoint:

    - Identifies common samples across all omics datasets
    - Reports omics-specific samples (measured in only some layers)
    - Creates a sample mapping file for downstream integration
    - Provides statistics on sample overlap

    This is a critical QC step before running MOFA2 or DIABLO, as both
    methods require matched samples. The output helps decide whether to:
    - Use only common samples (complete cases)
    - Use MOFA2's missing value handling (allows partial data)
    - Subset specific omics combinations
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Multi-Omics Sample Matching",
            workflow_type="multiomics_sample_matching",
            sample_id=request.sample_id,
            input_files={"omics_tables": request.omics_tables},
            parameters={},
        ),
    )

    background_tasks.add_task(
        multiomics_integrator.match_samples_across_omics,
        workflow.id,
        request.omics_tables,
        request.output_file,
        db,
    )

    return MultiOmicsResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"Sample matching across {len(request.omics_tables)} omics layers queued",
    )


@router.post("/complete-pipeline", response_model=MultiOmicsResponse)
async def run_complete_integration_pipeline(
    request: MultiOmicsWorkflowRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run complete multi-omics integration pipeline.

    This comprehensive workflow orchestrates the full multi-omics analysis:

    1. **Sample matching**: Identify common samples across omics
    2. **Data integration**: Run MOFA2 (unsupervised) or DIABLO (supervised)
    3. **Pathway enrichment**: Functional interpretation of results

    **When to use MOFA2 vs DIABLO:**

    Use MOFA2 when:
    - Exploratory analysis without predefined groups
    - Discovering sources of variation
    - Handling missing data across omics
    - Focus on biological interpretation

    Use DIABLO when:
    - Have defined phenotype/class labels
    - Interested in biomarker discovery
    - Building predictive models
    - Maximizing group discrimination

    The pipeline enables comprehensive understanding of multi-omics
    data from molecular measurements to biological insights.
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"Complete Multi-Omics Pipeline ({request.integration_method})",
            workflow_type="multiomics_complete",
            sample_id=request.sample_id,
            input_files={
                "data_matrices": request.data_matrices,
                "phenotype": request.phenotype_file,
            },
            parameters={
                "integration_method": request.integration_method,
                "run_pathway_enrichment": request.run_pathway_enrichment,
                "organism": request.organism,
            },
        ),
    )

    # Note: A complete orchestration would require sequential execution
    # Full implementation would use workflow engine (Nextflow/CWL)

    return MultiOmicsResponse(
        workflow_id=workflow.id,
        status="created",
        message=f"Complete multi-omics pipeline ({request.integration_method}) created. Use workflow engine for execution.",
    )
