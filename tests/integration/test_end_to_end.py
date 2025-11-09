"""
Integration tests for upload → process → retrieve workflow.

These tests verify the complete data processing pipeline:
1. Upload raw file
2. Process to unified format
3. Retrieve processed data
4. Verify metadata linking
5. Test cascade delete
"""

import io
import json
import os
import pytest
import requests
from pathlib import Path

# Configuration
BASE_URL = os.getenv("API_BASE_URL", "http://localhost:8001")
API_URL = f"{BASE_URL}/api/v1"

# Test credentials (should be environment variables in production)
import time

TEST_EMAIL = f"test_{int(time.time())}@example.com"  # Unique email per run
TEST_PASSWORD = "testpassword123"

# Test data directory
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"


class TestClient:
    """Test client with authentication."""

    def __init__(self):
        self.session = requests.Session()
        self.token = None
        self.user_id = None

    def register_and_login(self):
        """Register and login test user."""
        # Try to register (may fail if user exists)
        register_data = {
            "email": TEST_EMAIL,
            "password": TEST_PASSWORD,
            "full_name": "Test User",
        }

        response = self.session.post(f"{API_URL}/auth/register", json=register_data)

        if response.status_code in (200, 201):
            user_data = response.json()
            self.user_id = user_data.get("id")
            print(f"✓ Test user registered (id={self.user_id})")
        elif response.status_code == 400:
            # User already exists
            print(f"ℹ User already registered")
        else:
            print(f"⚠ Registration returned {response.status_code}: {response.text}")

        # Login
        login_data = {"username": TEST_EMAIL, "password": TEST_PASSWORD}
        response = self.session.post(
            f"{API_URL}/auth/login/access-token", data=login_data
        )

        if response.status_code != 200:
            print(f"Login failed with {response.status_code}: {response.text}")
            print(f"Attempted login with email: {TEST_EMAIL}")

        assert response.status_code == 200, f"Login failed: {response.text}"

        token_data = response.json()
        self.token = token_data["access_token"]
        self.session.headers.update({"Authorization": f"Bearer {self.token}"})

        print(f"✓ Logged in successfully")

    def create_project(self, name="Test Project", description="Integration test"):
        """Create a test project."""
        project_data = {"name": name, "description": description}
        response = self.session.post(f"{API_URL}/projects", json=project_data)
        assert response.status_code == 201, f"Project creation failed: {response.text}"
        project = response.json()
        print(f"✓ Created project {project['id']}: {project['name']}")
        return project

    def upload_file(
        self,
        project_id,
        file_path,
        process_file=True,
        sample_id=None,
        organism=None,
        reference_genome=None,
    ):
        """Upload a file with optional processing."""
        with open(file_path, "rb") as f:
            files = {"file": (file_path.name, f)}
            data = {"project_id": project_id, "process_file": str(process_file).lower()}
            if sample_id:
                data["sample_id"] = sample_id
            if organism:
                data["organism"] = organism
            if reference_genome:
                data["reference_genome"] = reference_genome

            response = self.session.post(
                f"{API_URL}/data/upload", files=files, data=data
            )

        assert response.status_code == 201, f"Upload failed: {response.text}"
        result = response.json()
        print(f"✓ Uploaded {file_path.name} → file_id={result['id']}")
        return result

    def get_processed_data(self, file_id):
        """Retrieve processed unified data."""
        response = self.session.get(f"{API_URL}/data/{file_id}/processed")
        return response

    def delete_project(self, project_id):
        """Delete project (cascade delete)."""
        response = self.session.delete(f"{API_URL}/projects/{project_id}")
        assert response.status_code == 204, f"Delete failed: {response.text}"
        print(f"✓ Deleted project {project_id}")


@pytest.fixture(scope="module")
def client():
    """Create authenticated test client."""
    c = TestClient()
    c.register_and_login()
    return c


@pytest.fixture
def test_project(client):
    """Create a test project for each test."""
    project = client.create_project()
    yield project
    # Cleanup
    try:
        client.delete_project(project["id"])
    except Exception as e:
        print(f"Cleanup warning: {e}")


class TestUploadProcessRetrieve:
    """Test upload → process → retrieve workflow."""

    def test_health_check(self):
        """Verify API is accessible."""
        response = requests.get(f"{BASE_URL}/healthz")
        assert response.status_code == 200
        print("✓ API health check passed")

    def test_transcriptomics_upload_and_process(self, client, test_project):
        """Test transcriptomics CSV upload and processing."""
        test_file = TEST_DATA_DIR / "sample_transcriptomics.csv"
        assert test_file.exists(), f"Test file not found: {test_file}"

        # Upload with processing
        result = client.upload_file(
            project_id=test_project["id"],
            file_path=test_file,
            process_file=True,
            sample_id="SAMPLE_001",
            organism="Homo sapiens",
        )

        # Verify upload result
        assert "id" in result

        # Check for processing info (may be "processing" or "processing_info")
        processing_info = result.get("processing") or result.get("processing_info", {})
        if not processing_info:
            print(f"⚠ Upload result keys: {result.keys()}")
            print(f"⚠ Full result: {result}")
        assert processing_info, "No processing info in upload result"

        # Should be processed
        assert processing_info.get("processed") == True
        assert processing_info.get("converted") == True
        assert "processed_file_id" in processing_info
        assert processing_info["omics_type"] == "transcriptomics"

        # Retrieve processed data (pass raw file ID to /processed endpoint)
        processed_file_id = processing_info.get("processed_file_id")
        assert processed_file_id, "No processed_file_id in processing info"

        raw_file_id = result["id"]
        response = client.get_processed_data(raw_file_id)
        assert response.status_code == 200

        processed_data = response.json()
        assert "unified_data" in processed_data
        unified = processed_data["unified_data"]

        # Verify unified format structure
        assert "metadata" in unified
        assert "records" in unified
        assert len(unified["records"]) > 0
        assert unified["metadata"]["omics_type"] == "transcriptomics"

        # Verify sample ID propagation
        assert unified["metadata"]["sample_id"] == "SAMPLE_001"
        assert unified["metadata"]["organism"] == "Homo sapiens"

        print(
            f"✓ Transcriptomics processing verified: {len(unified['records'])} records"
        )

    def test_genomics_upload_and_process(self, client, test_project):
        """Test genomics VCF upload and processing."""
        test_file = TEST_DATA_DIR / "sample_genomics.vcf"
        assert test_file.exists(), f"Test file not found: {test_file}"

        result = client.upload_file(
            project_id=test_project["id"],
            file_path=test_file,
            process_file=True,
            sample_id="SAMPLE_002",
            organism="Homo sapiens",
            reference_genome="hg38",
        )

        # Verify processing
        processing_info = result.get("processing") or result.get("processing_info", {})
        assert processing_info.get("processed") == True
        assert processing_info.get("converted") == True
        assert processing_info["omics_type"] == "genomics"

        # Retrieve processed data
        processed_file_id = processing_info.get("processed_file_id")
        assert processed_file_id, "No processed_file_id in processing info"

        raw_file_id = result["id"]
        response = client.get_processed_data(raw_file_id)
        assert response.status_code == 200

        unified = response.json()["unified_data"]
        assert unified["metadata"]["omics_type"] == "genomics"
        assert unified["metadata"]["reference_genome"] == "hg38"
        assert len(unified["records"]) >= 3  # At least 3 variants from test file

        print(f"✓ Genomics processing verified: {len(unified['records'])} variants")

    def test_proteomics_upload_and_process(self, client, test_project):
        """Test proteomics CSV upload and processing."""
        test_file = TEST_DATA_DIR / "sample_proteomics.csv"
        assert test_file.exists(), f"Test file not found: {test_file}"

        result = client.upload_file(
            project_id=test_project["id"],
            file_path=test_file,
            process_file=True,
            sample_id="SAMPLE_003",
        )

        # Verify processing
        processing_info = result.get("processing") or result.get("processing_info", {})
        assert processing_info.get("processed") == True
        assert processing_info.get("converted") == True
        assert processing_info["omics_type"] == "proteomics"

        # Retrieve processed data
        processed_file_id = processing_info.get("processed_file_id")
        assert processed_file_id, "No processed_file_id in processing info"
        raw_file_id = result["id"]
        response = client.get_processed_data(raw_file_id)
        assert response.status_code == 200

        unified = response.json()["unified_data"]
        assert unified["metadata"]["omics_type"] == "proteomics"
        assert len(unified["records"]) >= 4  # At least 4 proteins from test file

        print(f"✓ Proteomics processing verified: {len(unified['records'])} proteins")

    def test_upload_without_processing(self, client, test_project):
        """Test raw upload without processing."""
        test_file = TEST_DATA_DIR / "sample_transcriptomics.csv"

        result = client.upload_file(
            project_id=test_project["id"], file_path=test_file, process_file=False
        )

        # Should not be processed
        processing_info = result.get("processing") or result.get("processing_info", {})
        assert processing_info.get("processed") == False
        assert "processed_file_id" not in processing_info

        print("✓ Raw upload (no processing) verified")

    def test_cascade_delete(self, client):
        """Test cascade delete of project and associated files."""
        # Create project
        project = client.create_project(name="Delete Test Project")
        project_id = project["id"]

        # Upload multiple files
        test_file = TEST_DATA_DIR / "sample_transcriptomics.csv"
        for i in range(3):
            client.upload_file(
                project_id=project_id,
                file_path=test_file,
                process_file=True,
                sample_id=f"SAMPLE_{i}",
            )

        # Delete project (should cascade delete all files)
        client.delete_project(project_id)

        # Verify project is deleted
        response = client.session.get(f"{API_URL}/projects/{project_id}")
        assert response.status_code == 404

        print("✓ Cascade delete verified")

    def test_metadata_linking(self, client, test_project):
        """Test that raw and processed files are properly linked."""
        test_file = TEST_DATA_DIR / "sample_transcriptomics.csv"

        result = client.upload_file(
            project_id=test_project["id"], file_path=test_file, process_file=True
        )

        raw_file_id = result["id"]
        processing_info = result.get("processing") or result.get("processing_info", {})
        processed_file_id = processing_info.get("processed_file_id")

        assert processed_file_id, "No processed_file_id in processing info"

        # Verify files are different
        assert raw_file_id != processed_file_id

        # Get processed data via raw file ID (endpoint uses raw ID to find processed)
        response = client.get_processed_data(raw_file_id)
        if response.status_code != 200:
            print(f"❌ Get processed failed: {response.status_code} - {response.text}")
        assert response.status_code == 200

        print(
            f"✓ Metadata linking verified: raw={raw_file_id} → processed={processed_file_id}"
        )


def run_integration_tests():
    """Run integration tests with pytest."""
    import sys

    # Check if services are running
    try:
        response = requests.get(f"{BASE_URL}/healthz", timeout=5)
        if response.status_code != 200:
            print(f"❌ Backend health check failed: {response.status_code}")
            sys.exit(1)
    except Exception as e:
        print(f"❌ Cannot connect to backend at {BASE_URL}: {e}")
        print("Please ensure services are running with: ./docker-start.sh")
        sys.exit(1)

    print(f"\n{'='*60}")
    print("Running Integration Tests")
    print(f"API Base URL: {BASE_URL}")
    print(f"{'='*60}\n")

    # Run pytest
    exit_code = pytest.main(
        [__file__, "-v", "--tb=short", "-s"]  # Show print statements
    )

    return exit_code


if __name__ == "__main__":
    exit(run_integration_tests())
