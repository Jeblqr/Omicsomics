"""Tests for multi-omics integration API endpoints."""

import pytest
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_mofa2_integration(async_client: AsyncClient, auth_headers: dict):
    """Test MOFA2 integration endpoint."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
            "metabolomics": "/data/metabolite_features.csv",
        },
        "output_dir": "/output/mofa2",
        "n_factors": 10,
        "convergence_mode": "fast",
    }

    response = await async_client.post(
        "/api/multiomics/mofa2", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] in ["queued", "running", "completed"]
    assert "3 omics layers" in data["message"]


@pytest.mark.asyncio
async def test_diablo_integration(async_client: AsyncClient, auth_headers: dict):
    """Test DIABLO integration endpoint."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
        },
        "phenotype_file": "/data/phenotype.csv",
        "output_dir": "/output/diablo",
        "n_components": 2,
        "design_correlation": 0.1,
    }

    response = await async_client.post(
        "/api/multiomics/diablo", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] in ["queued", "running", "completed"]


@pytest.mark.asyncio
async def test_pathway_enrichment(async_client: AsyncClient, auth_headers: dict):
    """Test multi-omics pathway enrichment endpoint."""
    request_data = {
        "sample_id": 1,
        "feature_lists": {
            "transcriptomics": "/data/deg_list.csv",
            "proteomics": "/data/dap_list.csv",
        },
        "organism": "hsapiens",
        "output_dir": "/output/pathways",
    }

    response = await async_client.post(
        "/api/multiomics/pathway-enrichment", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] in ["queued", "running", "completed"]


@pytest.mark.asyncio
async def test_sample_matching(async_client: AsyncClient, auth_headers: dict):
    """Test sample matching endpoint."""
    request_data = {
        "sample_id": 1,
        "omics_tables": {
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
            "metabolomics": "/data/metabolite_features.csv",
        },
        "output_file": "/output/sample_mapping.json",
    }

    response = await async_client.post(
        "/api/multiomics/match-samples", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] in ["queued", "running", "completed"]


@pytest.mark.asyncio
async def test_complete_integration_pipeline_mofa2(
    async_client: AsyncClient, auth_headers: dict
):
    """Test complete multi-omics pipeline with MOFA2."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
        },
        "output_dir": "/output/complete_mofa2",
        "integration_method": "mofa2",
        "run_pathway_enrichment": True,
        "organism": "hsapiens",
    }

    response = await async_client.post(
        "/api/multiomics/complete-pipeline", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] in ["created", "queued"]
    assert "mofa2" in data["message"]


@pytest.mark.asyncio
async def test_complete_integration_pipeline_diablo(
    async_client: AsyncClient, auth_headers: dict
):
    """Test complete multi-omics pipeline with DIABLO."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
        },
        "phenotype_file": "/data/phenotype.csv",
        "output_dir": "/output/complete_diablo",
        "integration_method": "diablo",
        "run_pathway_enrichment": True,
        "organism": "hsapiens",
    }

    response = await async_client.post(
        "/api/multiomics/complete-pipeline", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "diablo" in data["message"]


@pytest.mark.asyncio
async def test_mofa2_with_custom_parameters(
    async_client: AsyncClient, auth_headers: dict
):
    """Test MOFA2 with custom parameters."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "proteomics": "/data/protein_abundance.csv",
        },
        "output_dir": "/output/mofa2_custom",
        "n_factors": 15,
        "convergence_mode": "slow",
    }

    response = await async_client.post(
        "/api/multiomics/mofa2", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data


@pytest.mark.asyncio
async def test_diablo_with_custom_design(async_client: AsyncClient, auth_headers: dict):
    """Test DIABLO with custom design correlation."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "metabolomics": "/data/metabolite_features.csv",
        },
        "phenotype_file": "/data/phenotype.csv",
        "output_dir": "/output/diablo_custom",
        "n_components": 3,
        "design_correlation": 0.5,
    }

    response = await async_client.post(
        "/api/multiomics/diablo", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data


@pytest.mark.asyncio
async def test_two_omics_integration(async_client: AsyncClient, auth_headers: dict):
    """Test integration with two omics layers."""
    request_data = {
        "sample_id": 1,
        "data_matrices": {
            "transcriptomics": "/data/rna_counts.csv",
            "metabolomics": "/data/metabolite_features.csv",
        },
        "output_dir": "/output/two_omics",
        "n_factors": 5,
    }

    response = await async_client.post(
        "/api/multiomics/mofa2", json=request_data, headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "2 omics layers" in data["message"]
