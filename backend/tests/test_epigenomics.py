"""Test epigenomics endpoints."""
from __future__ import annotations

import pytest
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_epigenomics_alignment_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test epigenomics alignment endpoint."""
    response = await async_client.post(
        "/api/v1/epigenomics/align",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "input_files": ["/tmp/test_R1.fastq.gz", "/tmp/test_R2.fastq.gz"],
            "reference_genome": "/tmp/genome/hg38",
            "output_bam": "/tmp/output/aligned.bam",
            "aligner": "bowtie2",
            "threads": 8,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "alignment" in data["message"].lower()


@pytest.mark.asyncio
async def test_epigenomics_alignment_unauthorized(async_client: AsyncClient):
    """Test epigenomics alignment without auth fails."""
    response = await async_client.post(
        "/api/v1/epigenomics/align",
        json={
            "sample_id": 1,
            "input_files": ["/tmp/test_R1.fastq.gz"],
            "reference_genome": "/tmp/genome/hg38",
            "output_bam": "/tmp/output/aligned.bam",
        },
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_peak_calling_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test peak calling endpoint."""
    response = await async_client.post(
        "/api/v1/epigenomics/peak-calling",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "treatment_bam": "/tmp/treatment.bam",
            "control_bam": "/tmp/control.bam",
            "output_dir": "/tmp/peaks",
            "peak_caller": "macs2",
            "genome_size": "hs",
            "broad": False,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "peak" in data["message"].lower()


@pytest.mark.asyncio
async def test_peak_calling_without_control(async_client: AsyncClient, auth_headers, test_project):
    """Test peak calling without control sample."""
    response = await async_client.post(
        "/api/v1/epigenomics/peak-calling",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "treatment_bam": "/tmp/treatment.bam",
            "output_dir": "/tmp/peaks",
            "peak_caller": "macs2",
            "genome_size": "mm",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"


@pytest.mark.asyncio
async def test_motif_analysis_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test motif analysis endpoint."""
    response = await async_client.post(
        "/api/v1/epigenomics/motif-analysis",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "peak_file": "/tmp/peaks.narrowPeak",
            "genome_fasta": "/tmp/genome/hg38.fa",
            "output_dir": "/tmp/motifs",
            "tool": "homer",
            "motif_length": [8, 10, 12],
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "motif" in data["message"].lower()


@pytest.mark.asyncio
async def test_bigwig_generation_endpoint(async_client: AsyncClient, auth_headers, test_project):
    """Test BigWig generation endpoint."""
    response = await async_client.post(
        "/api/v1/epigenomics/bigwig",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "input_bam": "/tmp/aligned.bam",
            "output_bigwig": "/tmp/coverage.bw",
            "genome_sizes": "/tmp/genome/hg38.chrom.sizes",
            "normalize": True,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"
    assert "bigwig" in data["message"].lower()


@pytest.mark.asyncio
async def test_complete_epigenomics_pipeline(async_client: AsyncClient, auth_headers, test_project):
    """Test complete epigenomics pipeline endpoint."""
    response = await async_client.post(
        "/api/v1/epigenomics/complete-pipeline",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "treatment_fastq": ["/tmp/treatment_R1.fastq.gz", "/tmp/treatment_R2.fastq.gz"],
            "control_fastq": ["/tmp/control_R1.fastq.gz", "/tmp/control_R2.fastq.gz"],
            "reference_genome": "/tmp/genome/hg38",
            "genome_fasta": "/tmp/genome/hg38.fa",
            "genome_sizes": "/tmp/genome/hg38.chrom.sizes",
            "output_dir": "/tmp/epigenomics_results",
            "aligner": "bowtie2",
            "peak_caller": "macs2",
            "genome_size": "hs",
            "run_motif_analysis": True,
            "generate_bigwig": True,
            "threads": 8,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "created"
    assert "complete" in data["message"].lower()


@pytest.mark.asyncio
async def test_epigenomics_with_broad_peaks(async_client: AsyncClient, auth_headers, test_project):
    """Test peak calling with broad peak mode."""
    response = await async_client.post(
        "/api/v1/epigenomics/peak-calling",
        headers=auth_headers,
        json={
            "sample_id": 1,
            "treatment_bam": "/tmp/treatment.bam",
            "control_bam": "/tmp/control.bam",
            "output_dir": "/tmp/peaks",
            "peak_caller": "macs2",
            "genome_size": "hs",
            "broad": True,
            "q_value": 0.05,
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert "workflow_id" in data
    assert data["status"] == "queued"


@pytest.mark.asyncio
async def test_epigenomics_invalid_sample(async_client: AsyncClient, auth_headers):
    """Test epigenomics endpoint with invalid sample ID."""
    response = await async_client.post(
        "/api/v1/epigenomics/align",
        headers=auth_headers,
        json={
            "sample_id": 99999,
            "input_files": ["/tmp/test.fastq.gz"],
            "reference_genome": "/tmp/genome/hg38",
            "output_bam": "/tmp/output/aligned.bam",
        },
    )
    # Should still create workflow record, but may fail during execution
    # The API accepts the request, validation happens later
    assert response.status_code == 200
