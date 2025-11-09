"""
Integration tests for archive processing.

Tests the complete workflow:
1. Upload archive (ZIP/TAR.GZ)
2. Preview contents
3. Select and extract file
4. Process to unified format
"""

import io
import zipfile
import tarfile
from pathlib import Path

import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession


@pytest.fixture
def sample_csv_content():
    """Sample CSV proteomics data."""
    return b"""protein_id,abundance,peptides
P12345,1500.5,3
Q67890,2300.8,5
R24680,890.2,2
"""


@pytest.fixture
def sample_vcf_content():
    """Sample VCF genomics data."""
    return b"""##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	12345	rs123	A	G	100	PASS	AF=0.5
chr1	67890	rs456	C	T	99	PASS	AF=0.3
"""


@pytest.fixture
def sample_json_content():
    """Sample JSON transcriptomics data."""
    return b"""{
  "sample_id": "test_sample",
  "genes": [
    {"gene_id": "ENSG001", "expression": 45.2, "tpm": 12.5},
    {"gene_id": "ENSG002", "expression": 102.8, "tpm": 28.3}
  ]
}
"""


@pytest.fixture
def test_zip_archive(sample_csv_content, sample_vcf_content, sample_json_content):
    """Create test ZIP archive with multiple omics files."""
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("proteomics/proteins.csv", sample_csv_content)
        zf.writestr("genomics/variants.vcf", sample_vcf_content)
        zf.writestr("transcriptomics/data.json", sample_json_content)
        zf.writestr("README.txt", b"Sample omics data archive")
    
    zip_buffer.seek(0)
    return zip_buffer


@pytest.fixture
def test_tar_gz_archive(sample_csv_content, sample_vcf_content):
    """Create test TAR.GZ archive."""
    import tempfile
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        # Create temp files
        csv_file = tmp_path / "proteins.csv"
        csv_file.write_bytes(sample_csv_content)
        
        vcf_file = tmp_path / "variants.vcf"
        vcf_file.write_bytes(sample_vcf_content)
        
        # Create TAR.GZ
        tar_buffer = io.BytesIO()
        with tarfile.open(fileobj=tar_buffer, mode='w:gz') as tf:
            tf.add(csv_file, arcname="proteomics/proteins.csv")
            tf.add(vcf_file, arcname="genomics/variants.vcf")
        
        tar_buffer.seek(0)
        return tar_buffer


class TestArchivePreview:
    """Test archive preview functionality."""
    
    @pytest.mark.asyncio
    async def test_preview_zip_archive(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_zip_archive: io.BytesIO,
    ):
        """Test previewing ZIP archive contents."""
        # Upload archive
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",  # Don't process the archive itself
            },
            files={"file": ("test_archive.zip", test_zip_archive, "application/zip")},
        )
        assert upload_response.status_code == 201
        datafile_id = upload_response.json()["id"]
        
        # Preview archive
        preview_response = await async_client.get(
            f"/api/v1/data/{datafile_id}/preview",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert preview_response.status_code == 200
        
        preview_data = preview_response.json()
        assert preview_data["is_archive"] is True
        assert preview_data["total_files"] == 4  # 3 data files + 1 README
        assert preview_data["archive_filename"] == "test_archive.zip"
        
        # Check files
        files = preview_data["files"]
        file_paths = {f["path"] for f in files}
        assert "proteomics/proteins.csv" in file_paths
        assert "genomics/variants.vcf" in file_paths
        assert "transcriptomics/data.json" in file_paths
        assert "README.txt" in file_paths
        
        # Check file metadata
        csv_file = next(f for f in files if f["name"] == "proteins.csv")
        assert csv_file["mime_type"] == "text/csv"
        assert csv_file["size"] > 0
    
    @pytest.mark.asyncio
    async def test_preview_tar_gz_archive(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_tar_gz_archive: io.BytesIO,
    ):
        """Test previewing TAR.GZ archive contents."""
        # Upload archive
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("test_archive.tar.gz", test_tar_gz_archive, "application/gzip")},
        )
        assert upload_response.status_code == 201
        datafile_id = upload_response.json()["id"]
        
        # Preview archive
        preview_response = await async_client.get(
            f"/api/v1/data/{datafile_id}/preview",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert preview_response.status_code == 200
        
        preview_data = preview_response.json()
        assert preview_data["is_archive"] is True
        assert preview_data["total_files"] == 2
        
        file_paths = {f["path"] for f in preview_data["files"]}
        assert "proteomics/proteins.csv" in file_paths
        assert "genomics/variants.vcf" in file_paths
    
    @pytest.mark.asyncio
    async def test_preview_non_archive_file(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        sample_csv_content: bytes,
    ):
        """Test previewing non-archive file returns appropriate message."""
        # Upload regular CSV
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("proteins.csv", io.BytesIO(sample_csv_content), "text/csv")},
        )
        assert upload_response.status_code == 201
        datafile_id = upload_response.json()["id"]
        
        # Try to preview
        preview_response = await async_client.get(
            f"/api/v1/data/{datafile_id}/preview",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert preview_response.status_code == 200
        
        preview_data = preview_response.json()
        assert preview_data["is_archive"] is False
        assert "message" in preview_data


class TestArchiveExtraction:
    """Test extracting and processing files from archives."""
    
    @pytest.mark.asyncio
    async def test_extract_and_process_csv_from_zip(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_zip_archive: io.BytesIO,
    ):
        """Test extracting and processing CSV file from ZIP."""
        # Upload archive
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("omics_data.zip", test_zip_archive, "application/zip")},
        )
        assert upload_response.status_code == 201
        archive_id = upload_response.json()["id"]
        
        # Extract and process CSV
        process_response = await async_client.post(
            f"/api/v1/data/{archive_id}/process-from-archive",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={"file_path": "proteomics/proteins.csv"},
        )
        assert process_response.status_code == 201
        
        result = process_response.json()
        assert result["extracted_filename"] == "proteins.csv"
        assert "datafile_id" in result
        assert "processed_file_id" in result
        assert result["archive_source"]["archive_id"] == archive_id
        
        # Verify processing info
        processing_info = result["processing_info"]
        assert processing_info["success"] is True
        assert "format_detected" in processing_info
        assert "features_count" in processing_info
        
        # Get processed data
        processed_response = await async_client.get(
            f"/api/v1/data/{result['datafile_id']}/processed",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert processed_response.status_code == 200
        
        unified_data = processed_response.json()
        assert "unified_data" in unified_data
        assert unified_data["unified_data"]["data_type"] == "proteomics"
    
    @pytest.mark.asyncio
    async def test_extract_and_process_vcf_from_zip(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_zip_archive: io.BytesIO,
    ):
        """Test extracting and processing VCF file from ZIP."""
        # Upload archive
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("omics_data.zip", test_zip_archive, "application/zip")},
        )
        assert upload_response.status_code == 201
        archive_id = upload_response.json()["id"]
        
        # Extract and process VCF
        process_response = await async_client.post(
            f"/api/v1/data/{archive_id}/process-from-archive",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={"file_path": "genomics/variants.vcf"},
        )
        assert process_response.status_code == 201
        
        result = process_response.json()
        assert result["extracted_filename"] == "variants.vcf"
        
        # Get processed data
        processed_response = await async_client.get(
            f"/api/v1/data/{result['datafile_id']}/processed",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert processed_response.status_code == 200
        
        unified_data = processed_response.json()
        assert unified_data["unified_data"]["data_type"] == "genomics"
    
    @pytest.mark.asyncio
    async def test_extract_nonexistent_file_fails(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_zip_archive: io.BytesIO,
    ):
        """Test extracting non-existent file returns error."""
        # Upload archive
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("omics_data.zip", test_zip_archive, "application/zip")},
        )
        assert upload_response.status_code == 201
        archive_id = upload_response.json()["id"]
        
        # Try to extract non-existent file
        process_response = await async_client.post(
            f"/api/v1/data/{archive_id}/process-from-archive",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={"file_path": "nonexistent/file.txt"},
        )
        assert process_response.status_code in [400, 404, 500]


class TestArchiveWorkflow:
    """Test complete archive processing workflow."""
    
    @pytest.mark.asyncio
    async def test_complete_workflow_zip(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_zip_archive: io.BytesIO,
    ):
        """Test complete workflow: upload → preview → extract → process."""
        # Step 1: Upload archive
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("omics_data.zip", test_zip_archive, "application/zip")},
        )
        assert upload_response.status_code == 201
        archive_id = upload_response.json()["id"]
        
        # Step 2: Preview contents
        preview_response = await async_client.get(
            f"/api/v1/data/{archive_id}/preview",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert preview_response.status_code == 200
        preview_data = preview_response.json()
        assert len(preview_data["files"]) > 0
        
        # Step 3: Process multiple files
        processed_files = []
        
        for file_info in preview_data["files"]:
            # Skip README
            if file_info["name"] == "README.txt":
                continue
            
            # Extract and process
            process_response = await async_client.post(
                f"/api/v1/data/{archive_id}/process-from-archive",
                headers={"Authorization": f"Bearer {test_user_token}"},
                data={"file_path": file_info["path"]},
            )
            
            if process_response.status_code == 201:
                result = process_response.json()
                processed_files.append(result)
        
        # Verify all data files were processed
        assert len(processed_files) >= 2  # At least CSV and VCF
        
        # Verify each has processed data
        for processed_file in processed_files:
            get_response = await async_client.get(
                f"/api/v1/data/{processed_file['datafile_id']}/processed",
                headers={"Authorization": f"Bearer {test_user_token}"},
            )
            assert get_response.status_code == 200
            
            unified_data = get_response.json()
            assert "unified_data" in unified_data
            assert unified_data["unified_data"]["data_type"] in [
                "proteomics",
                "genomics",
                "transcriptomics",
            ]
    
    @pytest.mark.asyncio
    async def test_complete_workflow_tar_gz(
        self,
        async_client: AsyncClient,
        test_user_token: str,
        test_project: dict,
        test_tar_gz_archive: io.BytesIO,
    ):
        """Test complete workflow with TAR.GZ archive."""
        # Upload
        upload_response = await async_client.post(
            "/api/v1/data/upload",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={
                "project_id": test_project["id"],
                "process_file": "false",
            },
            files={"file": ("omics_data.tar.gz", test_tar_gz_archive, "application/gzip")},
        )
        assert upload_response.status_code == 201
        archive_id = upload_response.json()["id"]
        
        # Preview
        preview_response = await async_client.get(
            f"/api/v1/data/{archive_id}/preview",
            headers={"Authorization": f"Bearer {test_user_token}"},
        )
        assert preview_response.status_code == 200
        
        # Process first file
        preview_data = preview_response.json()
        first_file = preview_data["files"][0]
        
        process_response = await async_client.post(
            f"/api/v1/data/{archive_id}/process-from-archive",
            headers={"Authorization": f"Bearer {test_user_token}"},
            data={"file_path": first_file["path"]},
        )
        assert process_response.status_code == 201
