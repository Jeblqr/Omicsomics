"""Test proteomics endpoints."""
from __future__ import annotations

import pytest
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_raw_conversion_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test raw to mzML conversion endpoint."""
    response = await async_client.post(
        "/api/v1/proteomics/convert-raw",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "raw_files": ["/tmp/sample1.raw", "/tmp/sample2.raw"],
            "output_dir": "/tmp/proteomics/mzml",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "conversion" in data["message"].lower()


@pytest.mark.asyncio
async def test_raw_conversion_unauthorized(async_client: AsyncClient):
    """Test raw conversion without auth fails."""
    response = await async_client.post(
        "/api/v1/proteomics/convert-raw",
        json={
            "sample_id": 1,
            "raw_files": ["/tmp/sample.raw"],
            "output_dir": "/tmp/mzml",
        },
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_maxquant_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test MaxQuant analysis endpoint."""
    response = await async_client.post(
        "/api/v1/proteomics/maxquant",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "raw_files": ["/tmp/sample1.raw", "/tmp/sample2.raw"],
            "fasta_file": "/tmp/uniprot_human.fasta",
            "output_dir": "/tmp/proteomics/maxquant",
            "lfq": True,
            "match_between_runs": True,
            "max_missed_cleavages": 2,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "maxquant" in data["message"].lower()


@pytest.mark.asyncio
async def test_maxquant_without_lfq(async_client: AsyncClient, auth_headers, test_project):
    """Test MaxQuant without LFQ."""
    response = await async_client.post(
        "/api/v1/proteomics/maxquant",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "raw_files": ["/tmp/sample.raw"],
            "fasta_file": "/tmp/database.fasta",
            "output_dir": "/tmp/output",
            "lfq": False,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"


@pytest.mark.asyncio
async def test_msfragger_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test MSFragger analysis endpoint."""
    response = await async_client.post(
        "/api/v1/proteomics/msfragger",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample1.mzML", "/tmp/sample2.mzML"],
            "fasta_file": "/tmp/uniprot_human.fasta",
            "output_dir": "/tmp/proteomics/msfragger",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "msfragger" in data["message"].lower()


@pytest.mark.asyncio
async def test_lfq_quantification_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test LFQ quantification endpoint."""
    response = await async_client.post(
        "/api/v1/proteomics/lfq-quantification",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "protein_groups_file": "/tmp/maxquant/proteinGroups.txt",
            "output_file": "/tmp/proteomics/lfq_quantification.csv",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "quantification" in data["message"].lower()


@pytest.mark.asyncio
async def test_complete_proteomics_pipeline(async_client: AsyncClient, auth_headers, test_project):
    """Test complete proteomics pipeline endpoint."""
    response = await async_client.post(
        "/api/v1/proteomics/complete-pipeline",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "raw_files": ["/tmp/sample1.raw", "/tmp/sample2.raw"],
            "fasta_file": "/tmp/uniprot_human.fasta",
            "output_dir": "/tmp/proteomics_results",
            "tool": "maxquant",
            "lfq": True,
            "convert_raw": True,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "created"
    assert "complete" in data["message"].lower()


@pytest.mark.asyncio
async def test_complete_pipeline_with_msfragger(async_client: AsyncClient, auth_headers, test_project):
    """Test complete proteomics pipeline with MSFragger."""
    response = await async_client.post(
        "/api/v1/proteomics/complete-pipeline",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "raw_files": ["/tmp/sample.mzML"],
            "fasta_file": "/tmp/database.fasta",
            "output_dir": "/tmp/results",
            "tool": "msfragger",
            "lfq": False,
            "convert_raw": False,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "created"
    assert "msfragger" in data["message"].lower()


@pytest.mark.asyncio
async def test_proteomics_with_multiple_raw_files(async_client: AsyncClient, auth_headers, test_project):
    """Test proteomics with multiple raw files."""
    response = await async_client.post(
        "/api/v1/proteomics/maxquant",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "raw_files": [
                "/tmp/sample1.raw",
                "/tmp/sample2.raw",
                "/tmp/sample3.raw",
                "/tmp/sample4.raw",
            ],
            "fasta_file": "/tmp/database.fasta",
            "output_dir": "/tmp/output",
            "lfq": True,
            "match_between_runs": True,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data


@pytest.mark.asyncio
async def test_proteomics_invalid_sample(async_client: AsyncClient, auth_headers):
    """Test proteomics endpoint with invalid sample ID."""
    response = await async_client.post(
        "/api/v1/proteomics/maxquant",
        headers=auth_headers,
        json={
            "sample_id": 99999,
            "raw_files": ["/tmp/sample.raw"],
            "fasta_file": "/tmp/database.fasta",
            "output_dir": "/tmp/output",
        },
    )
    # Should still create workflow record, but may fail during execution
    assert response.status_code == 200
