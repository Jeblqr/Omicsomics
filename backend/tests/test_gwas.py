"""Tests for GWAS API endpoints."""

import pytest
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_plink_qc(async_client: AsyncClient, auth_headers: dict):
    """Test PLINK QC endpoint."""
    request_data = {
        "sample_id": 1,
        "bed_file": "/data/genotypes.bed",
        "bim_file": "/data/genotypes.bim",
        "fam_file": "/data/genotypes.fam",
        "output_prefix": "/output/qc_filtered",
        "geno": 0.02,
        "mind": 0.02,
        "maf": 0.01,
        "hwe": 1e-6,
    }

    response = await async_client.post("/api/gwas/qc", json=request_data, headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] in ["queued", "running", "completed"]


@pytest.mark.asyncio
async def test_association_test_quantitative(async_client: AsyncClient, auth_headers: dict):
    """Test GWAS association test for quantitative trait."""
    request_data = {
        "sample_id": 1,
        "bed_file": "/data/qc_filtered.bed",
        "phenotype_file": "/data/height.phen",
        "output_prefix": "/output/gwas_height",
        "binary_trait": False,
    }

    response = await async_client.post(
        "/api/gwas/association", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "linear" in data["message"]


@pytest.mark.asyncio
async def test_association_test_binary(async_client: AsyncClient, auth_headers: dict):
    """Test GWAS association test for binary trait."""
    request_data = {
        "sample_id": 1,
        "bed_file": "/data/qc_filtered.bed",
        "phenotype_file": "/data/disease.phen",
        "output_prefix": "/output/gwas_disease",
        "covariates_file": "/data/covariates.txt",
        "binary_trait": True,
    }

    response = await async_client.post(
        "/api/gwas/association", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "logistic" in data["message"]


@pytest.mark.asyncio
async def test_ld_calculation(async_client: AsyncClient, auth_headers: dict):
    """Test LD calculation endpoint."""
    request_data = {
        "sample_id": 1,
        "bed_file": "/data/genotypes.bed",
        "output_prefix": "/output/ld_matrix",
        "ld_window": 1000,
        "ld_window_r2": 0.2,
    }

    response = await async_client.post("/api/gwas/ld", json=request_data, headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data


@pytest.mark.asyncio
async def test_prs_calculation(async_client: AsyncClient, auth_headers: dict):
    """Test PRS calculation endpoint."""
    request_data = {
        "sample_id": 1,
        "bed_file": "/data/target_genotypes.bed",
        "weights_file": "/data/gwas_weights.txt",
        "output_prefix": "/output/prs_scores",
    }

    response = await async_client.post("/api/gwas/prs", json=request_data, headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "PRS" in data["message"]


@pytest.mark.asyncio
async def test_mtag_analysis(async_client: AsyncClient, auth_headers: dict):
    """Test MTAG cross-trait analysis endpoint."""
    request_data = {
        "sample_id": 1,
        "summary_stats_files": {
            "height": "/data/height_sumstats.txt",
            "bmi": "/data/bmi_sumstats.txt",
            "weight": "/data/weight_sumstats.txt",
        },
        "output_dir": "/output/mtag",
    }

    response = await async_client.post("/api/gwas/mtag", json=request_data, headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "3 traits" in data["message"]


@pytest.mark.asyncio
async def test_mtag_two_traits(async_client: AsyncClient, auth_headers: dict):
    """Test MTAG with two traits."""
    request_data = {
        "sample_id": 1,
        "summary_stats_files": {
            "trait1": "/data/trait1_sumstats.txt",
            "trait2": "/data/trait2_sumstats.txt",
        },
        "output_dir": "/output/mtag_two",
    }

    response = await async_client.post("/api/gwas/mtag", json=request_data, headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "2 traits" in data["message"]


@pytest.mark.asyncio
async def test_qc_with_custom_thresholds(async_client: AsyncClient, auth_headers: dict):
    """Test QC with custom thresholds."""
    request_data = {
        "sample_id": 1,
        "bed_file": "/data/genotypes.bed",
        "bim_file": "/data/genotypes.bim",
        "fam_file": "/data/genotypes.fam",
        "output_prefix": "/output/strict_qc",
        "geno": 0.01,
        "mind": 0.01,
        "maf": 0.05,
        "hwe": 1e-10,
    }

    response = await async_client.post("/api/gwas/qc", json=request_data, headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
