"""Test metabolomics endpoints."""
from __future__ import annotations

import pytest
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_feature_detection_xcms_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test XCMS feature detection endpoint."""
    response = await async_client.post(
        "/api/v1/metabolomics/feature-detection",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample1.mzML", "/tmp/sample2.mzML"],
            "output_dir": "/tmp/metabolomics/features",
            "tool": "xcms",
            "ppm": 25,
            "peakwidth_min": 10.0,
            "peakwidth_max": 60.0,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "xcms" in data["message"].lower()


@pytest.mark.asyncio
async def test_feature_detection_mzmine_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test MZmine feature detection endpoint."""
    response = await async_client.post(
        "/api/v1/metabolomics/feature-detection",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample1.mzML", "/tmp/sample2.mzML"],
            "output_dir": "/tmp/metabolomics/features",
            "tool": "mzmine",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "mzmine" in data["message"].lower()


@pytest.mark.asyncio
async def test_feature_detection_unauthorized(async_client: AsyncClient):
    """Test feature detection without auth fails."""
    response = await async_client.post(
        "/api/v1/metabolomics/feature-detection",
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample.mzML"],
            "output_dir": "/tmp/features",
            "tool": "xcms",
        },
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_feature_detection_invalid_tool(async_client: AsyncClient, auth_headers, test_project):
    """Test feature detection with invalid tool."""
    response = await async_client.post(
        "/api/v1/metabolomics/feature-detection",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample.mzML"],
            "output_dir": "/tmp/features",
            "tool": "invalid_tool",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "failed"
    assert "unsupported" in data["message"].lower()


@pytest.mark.asyncio
async def test_gnps_annotation_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test GNPS annotation endpoint."""
    response = await async_client.post(
        "/api/v1/metabolomics/spectral-annotation",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "input_file": "/tmp/spectra.mgf",
            "output_dir": "/tmp/metabolomics/annotations",
            "tool": "gnps",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "gnps" in data["message"].lower()


@pytest.mark.asyncio
async def test_msdial_annotation_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test MS-DIAL annotation endpoint."""
    response = await async_client.post(
        "/api/v1/metabolomics/spectral-annotation",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "input_file": "/tmp/features.csv",
            "output_dir": "/tmp/metabolomics/annotations",
            "tool": "msdial",
            "msp_library": "/tmp/library.msp",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "msdial" in data["message"].lower()


@pytest.mark.asyncio
async def test_msdial_annotation_without_library(async_client: AsyncClient, auth_headers, test_project):
    """Test MS-DIAL annotation without required library."""
    response = await async_client.post(
        "/api/v1/metabolomics/spectral-annotation",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "input_file": "/tmp/features.csv",
            "output_dir": "/tmp/annotations",
            "tool": "msdial",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "failed"
    assert "library" in data["message"].lower()


@pytest.mark.asyncio
async def test_quantification_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test feature quantification endpoint."""
    response = await async_client.post(
        "/api/v1/metabolomics/quantification",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "feature_table": "/tmp/features.csv",
            "sample_metadata": "/tmp/metadata.csv",
            "output_dir": "/tmp/metabolomics/quantification",
            "normalization": "median",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "quantification" in data["message"].lower()


@pytest.mark.asyncio
async def test_quantification_with_quantile_normalization(async_client: AsyncClient, auth_headers, test_project):
    """Test quantification with quantile normalization."""
    response = await async_client.post(
        "/api/v1/metabolomics/quantification",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "feature_table": "/tmp/features.csv",
            "sample_metadata": "/tmp/metadata.csv",
            "output_dir": "/tmp/quantification",
            "normalization": "quantile",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert "quantile" in data["message"].lower()


@pytest.mark.asyncio
async def test_complete_metabolomics_pipeline(async_client: AsyncClient, auth_headers, test_project):
    """Test complete metabolomics pipeline endpoint."""
    response = await async_client.post(
        "/api/v1/metabolomics/complete-pipeline",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample1.mzML", "/tmp/sample2.mzML", "/tmp/sample3.mzML"],
            "sample_metadata": "/tmp/metadata.csv",
            "output_dir": "/tmp/metabolomics_results",
            "feature_tool": "xcms",
            "annotation_tool": "gnps",
            "normalization": "median",
            "run_annotation": True,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "created"
    assert "complete" in data["message"].lower()


@pytest.mark.asyncio
async def test_complete_pipeline_with_mzmine(async_client: AsyncClient, auth_headers, test_project):
    """Test complete metabolomics pipeline with MZmine."""
    response = await async_client.post(
        "/api/v1/metabolomics/complete-pipeline",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": ["/tmp/sample.mzML"],
            "sample_metadata": "/tmp/metadata.csv",
            "output_dir": "/tmp/results",
            "feature_tool": "mzmine",
            "annotation_tool": "msdial",
            "normalization": "quantile",
            "run_annotation": False,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "created"
    assert "mzmine" in data["message"].lower()


@pytest.mark.asyncio
async def test_metabolomics_with_many_samples(async_client: AsyncClient, auth_headers, test_project):
    """Test metabolomics with multiple samples."""
    response = await async_client.post(
        "/api/v1/metabolomics/feature-detection",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "mzml_files": [
                "/tmp/sample1.mzML",
                "/tmp/sample2.mzML",
                "/tmp/sample3.mzML",
                "/tmp/sample4.mzML",
                "/tmp/sample5.mzML",
            ],
            "output_dir": "/tmp/features",
            "tool": "xcms",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data


@pytest.mark.asyncio
async def test_metabolomics_invalid_sample(async_client: AsyncClient, auth_headers):
    """Test metabolomics endpoint with invalid sample ID."""
    response = await async_client.post(
        "/api/v1/metabolomics/feature-detection",
        headers=auth_headers,
        json={
            "sample_id": 99999,
            "mzml_files": ["/tmp/sample.mzML"],
            "output_dir": "/tmp/features",
            "tool": "xcms",
        },
    )
    # Should still create workflow record, but may fail during execution
    assert response.status_code == 200
