"""GWAS (Genome-Wide Association Study) API endpoints."""

import logging
from typing import Annotated

from app.api.dependencies import get_current_user, get_db
from app.models.user import User
from app.pipelines.gwas import gwas_analyzer
from app.schemas import workflows as workflow_schema
from app.services import workflows as workflow_service
from fastapi import APIRouter, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from sqlalchemy.ext.asyncio import AsyncSession

logger = logging.getLogger(__name__)

router = APIRouter(tags=["gwas"])


# Request/Response Models
class PlinkQCRequest(BaseModel):
    """Request for PLINK quality control."""

    sample_id: int = Field(..., description="Sample ID")
    bed_file: str = Field(..., description="PLINK BED file path")
    bim_file: str = Field(..., description="PLINK BIM file path")
    fam_file: str = Field(..., description="PLINK FAM file path")
    output_prefix: str = Field(..., description="Output file prefix")
    geno: float = Field(default=0.02, description="SNP missing rate threshold")
    mind: float = Field(default=0.02, description="Individual missing rate threshold")
    maf: float = Field(default=0.01, description="Minor allele frequency threshold")
    hwe: float = Field(
        default=1e-6, description="Hardy-Weinberg equilibrium p-value threshold"
    )


class AssociationTestRequest(BaseModel):
    """Request for GWAS association test."""

    sample_id: int = Field(..., description="Sample ID")
    bed_file: str = Field(..., description="PLINK BED file (QC-filtered)")
    phenotype_file: str = Field(..., description="Phenotype file path")
    output_prefix: str = Field(..., description="Output file prefix")
    covariates_file: str | None = Field(
        default=None, description="Optional covariates file"
    )
    binary_trait: bool = Field(
        default=False, description="Is phenotype binary (case/control)?"
    )


class LDCalculationRequest(BaseModel):
    """Request for LD calculation."""

    sample_id: int = Field(..., description="Sample ID")
    bed_file: str = Field(..., description="PLINK BED file")
    output_prefix: str = Field(..., description="Output file prefix")
    ld_window: int = Field(default=1000, description="LD window size in kb")
    ld_window_r2: float = Field(default=0.2, description="LD r2 threshold")


class PRSRequest(BaseModel):
    """Request for Polygenic Risk Score calculation."""

    sample_id: int = Field(..., description="Sample ID")
    bed_file: str = Field(..., description="PLINK BED file (target genotypes)")
    weights_file: str = Field(..., description="SNP weights file (summary statistics)")
    output_prefix: str = Field(..., description="Output file prefix")


class MTAGRequest(BaseModel):
    """Request for Multi-Trait Analysis of GWAS."""

    sample_id: int = Field(..., description="Sample ID")
    summary_stats_files: dict[str, str] = Field(
        ...,
        description="Dictionary mapping trait names to summary statistics files",
        example={
            "height": "/data/height_sumstats.txt",
            "bmi": "/data/bmi_sumstats.txt",
        },
    )
    output_dir: str = Field(..., description="Output directory")


class GWASResponse(BaseModel):
    """Response for GWAS operations."""

    workflow_id: int
    status: str
    message: str


# Endpoints
@router.post("/qc", response_model=GWASResponse)
async def run_plink_qc(
    request: PlinkQCRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run PLINK quality control on genotype data.

    Quality control steps include:
    - Filtering SNPs with high missing rate (--geno)
    - Filtering individuals with high missing rate (--mind)
    - Filtering rare variants by MAF (--maf)
    - Hardy-Weinberg equilibrium test (--hwe)

    This is a critical preprocessing step before GWAS to ensure:
    - Data quality and reliability
    - Removal of genotyping errors
    - Compliance with HWE assumptions
    - Adequate sample and SNP coverage

    Output: QC-filtered PLINK binary files (.bed, .bim, .fam)
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="GWAS Quality Control (PLINK)",
            workflow_type="gwas_qc",
            sample_id=request.sample_id,
            input_files={
                "bed": request.bed_file,
                "bim": request.bim_file,
                "fam": request.fam_file,
            },
            parameters={
                "geno": request.geno,
                "mind": request.mind,
                "maf": request.maf,
                "hwe": request.hwe,
            },
        ),
    )

    params = {
        "geno": request.geno,
        "mind": request.mind,
        "maf": request.maf,
        "hwe": request.hwe,
    }

    background_tasks.add_task(
        gwas_analyzer.run_plink_qc,
        workflow.id,
        request.bed_file,
        request.bim_file,
        request.fam_file,
        request.output_prefix,
        params,
        db,
    )

    return GWASResponse(
        workflow_id=workflow.id,
        status="queued",
        message="PLINK QC queued for execution",
    )


@router.post("/association", response_model=GWASResponse)
async def run_association_test(
    request: AssociationTestRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run GWAS association test.

    Performs genome-wide association testing between genetic variants
    and phenotype using PLINK.

    Test types:
    - Linear regression: for quantitative traits
    - Logistic regression: for binary traits (case/control)

    Features:
    - Covariate adjustment (e.g., age, sex, PCs)
    - Multiple testing correction (Bonferroni, FDR)
    - Association statistics (beta, SE, p-value)

    Output: Association results with:
    - SNP identifier
    - Effect size (beta/OR)
    - Standard error
    - P-value
    - Adjusted p-values
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="GWAS Association Test",
            workflow_type="gwas_association",
            sample_id=request.sample_id,
            input_files={
                "bed": request.bed_file,
                "phenotype": request.phenotype_file,
                "covariates": request.covariates_file,
            },
            parameters={"binary_trait": request.binary_trait},
        ),
    )

    params = {"binary_trait": request.binary_trait}

    background_tasks.add_task(
        gwas_analyzer.run_association_test,
        workflow.id,
        request.bed_file,
        request.phenotype_file,
        request.output_prefix,
        request.covariates_file,
        params,
        db,
    )

    return GWASResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"GWAS {'logistic' if request.binary_trait else 'linear'} association test queued",
    )


@router.post("/ld", response_model=GWASResponse)
async def calculate_ld(
    request: LDCalculationRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Calculate linkage disequilibrium (LD) between SNPs.

    LD measures the non-random association of alleles at different loci.
    It's crucial for:
    - Fine-mapping causal variants
    - Understanding genetic architecture
    - LD pruning for PCA
    - Identifying independent signals

    Parameters:
    - LD window: Distance (kb) within which to calculate LD
    - r² threshold: Minimum r² to report

    Output: Pairwise LD matrix with:
    - SNP pairs
    - r² values (correlation)
    - D' values
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="LD Calculation",
            workflow_type="gwas_ld",
            sample_id=request.sample_id,
            input_files={"bed": request.bed_file},
            parameters={
                "ld_window": request.ld_window,
                "ld_window_r2": request.ld_window_r2,
            },
        ),
    )

    params = {
        "ld_window": request.ld_window,
        "ld_window_r2": request.ld_window_r2,
    }

    background_tasks.add_task(
        gwas_analyzer.calculate_ld,
        workflow.id,
        request.bed_file,
        request.output_prefix,
        params,
        db,
    )

    return GWASResponse(
        workflow_id=workflow.id,
        status="queued",
        message="LD calculation queued",
    )


@router.post("/prs", response_model=GWASResponse)
async def calculate_prs(
    request: PRSRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Calculate Polygenic Risk Score (PRS).

    PRS aggregates the effects of many genetic variants to predict
    an individual's genetic predisposition to a trait or disease.

    Applications:
    - Disease risk prediction
    - Stratified medicine
    - Drug response prediction
    - Genetic counseling

    Process:
    1. Use GWAS summary statistics as weights
    2. Calculate weighted sum of risk alleles per individual
    3. Normalize scores

    Output: PRS for each individual with:
    - Individual ID
    - Number of variants scored
    - Raw PRS
    - Standardized PRS
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name="Polygenic Risk Score Calculation",
            workflow_type="gwas_prs",
            sample_id=request.sample_id,
            input_files={
                "bed": request.bed_file,
                "weights": request.weights_file,
            },
            parameters={},
        ),
    )

    background_tasks.add_task(
        gwas_analyzer.run_prs_calculation,
        workflow.id,
        request.bed_file,
        request.weights_file,
        request.output_prefix,
        db,
    )

    return GWASResponse(
        workflow_id=workflow.id,
        status="queued",
        message="PRS calculation queued",
    )


@router.post("/mtag", response_model=GWASResponse)
async def run_mtag_analysis(
    request: MTAGRequest,
    background_tasks: BackgroundTasks,
    db: Annotated[AsyncSession, Depends(get_db)],
    current_user: Annotated[User, Depends(get_current_user)],
):
    """
    Run Multi-Trait Analysis of GWAS (MTAG).

    MTAG performs cross-trait meta-analysis to:
    - Boost power by leveraging genetic correlation
    - Identify shared genetic architecture
    - Improve effect size estimation
    - Discover pleiotropic variants

    When to use MTAG:
    - Multiple genetically correlated traits
    - Want to boost discovery power
    - Interest in shared biology

    Method:
    - Inverse-variance weighted meta-analysis
    - Accounts for sample overlap
    - Estimates genetic correlation
    - Performs cross-trait analysis

    Output:
    - Meta-analyzed summary statistics
    - Genetic correlation matrix
    - Improved p-values and effect sizes
    - Pleiotropic SNP identification
    """
    workflow = await workflow_service.create_workflow(
        db,
        workflow_schema.WorkflowCreate(
            name=f"MTAG Analysis ({len(request.summary_stats_files)} traits)",
            workflow_type="gwas_mtag",
            sample_id=request.sample_id,
            input_files={"summary_stats": request.summary_stats_files},
            parameters={},
        ),
    )

    background_tasks.add_task(
        gwas_analyzer.run_mtag_analysis,
        workflow.id,
        request.summary_stats_files,
        request.output_dir,
        None,
        db,
    )

    return GWASResponse(
        workflow_id=workflow.id,
        status="queued",
        message=f"MTAG analysis with {len(request.summary_stats_files)} traits queued",
    )
