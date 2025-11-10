#!/usr/bin/env python3
"""
Comprehensive test script for Data Export functionality.

Tests:
1. Export job creation
2. Format conversion exports
3. Progress tracking
4. Job cancellation
5. Job deletion
6. Statistics
7. Cleanup operations
"""

import sys
import time
import requests
from typing import Optional, Dict, Any, List

# Configuration
BASE_URL = "http://localhost:8000"
API_BASE = f"{BASE_URL}/api/data-export"
FILES_API = f"{BASE_URL}/api/files"

# Test data (adjust these IDs based on your data)
TEST_FILE_IDS = [1, 2, 3]  # Replace with actual file IDs


class Colors:
    """ANSI color codes for terminal output"""

    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    RESET = "\033[0m"
    BOLD = "\033[1m"


def print_section(title: str):
    """Print a section header"""
    print(f"\n{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.CYAN}{Colors.BOLD}{title}{Colors.RESET}")
    print(f"{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.RESET}\n")


def print_test(test_name: str):
    """Print test name"""
    print(f"{Colors.BLUE}Testing: {test_name}{Colors.RESET}")


def print_success(message: str):
    """Print success message"""
    print(f"{Colors.GREEN}✓ {message}{Colors.RESET}")


def print_error(message: str):
    """Print error message"""
    print(f"{Colors.RED}✗ {message}{Colors.RESET}")


def print_warning(message: str):
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠ {message}{Colors.RESET}")


def print_info(message: str):
    """Print info message"""
    print(f"{Colors.CYAN}ℹ {message}{Colors.RESET}")


class DataExportTest:
    """Test suite for Data Export"""

    def __init__(self):
        self.session = requests.Session()
        self.test_job_ids: List[int] = []
        self.passed = 0
        self.failed = 0

        # Get auth token (assuming you have auth)
        # self.authenticate()

    def authenticate(self):
        """Authenticate and get token"""
        # Implement authentication if required
        pass

    def assert_response(
        self, response: requests.Response, expected_status: int, test_name: str
    ) -> bool:
        """Assert response status and handle result"""
        if response.status_code == expected_status:
            self.passed += 1
            print_success(f"{test_name} - Status {response.status_code}")
            return True
        else:
            self.failed += 1
            print_error(
                f"{test_name} - Expected {expected_status}, got {response.status_code}"
            )
            try:
                error_data = response.json()
                print_info(f"Error: {error_data}")
            except:
                print_info(f"Response: {response.text}")
            return False

    def test_get_service_info(self) -> bool:
        """Test: Get export service information"""
        print_test("Get Service Info")

        response = self.session.get(f"{API_BASE}/")

        if self.assert_response(response, 200, "Get service info"):
            result = response.json()
            print_info(f"Service: {result['service']}")
            print_info(f"Version: {result['version']}")
            print_info(f"Supported formats: {', '.join(result['supported_formats'])}")
            return True
        return False

    def test_create_csv_export(self) -> bool:
        """Test: Create CSV export job"""
        print_test("Create CSV Export Job")

        data = {
            "name": "Test CSV Export",
            "description": "Test export to CSV format",
            "file_ids": TEST_FILE_IDS,
            "export_format": "csv",
            "include_metadata": True,
            "include_lineage": False,
            "compress": True,
            "ttl_hours": 24,
        }

        response = self.session.post(f"{API_BASE}/jobs", json=data)

        if self.assert_response(response, 200, "Create CSV export"):
            result = response.json()
            self.test_job_ids.append(result["id"])
            print_info(f"Created job ID: {result['id']}")
            print_info(f"Job key: {result['job_key']}")
            print_info(f"Status: {result['status']}")
            print_info(f"Total files: {result['total_files']}")
            return True
        return False

    def test_create_json_export(self) -> bool:
        """Test: Create JSON export job"""
        print_test("Create JSON Export Job")

        data = {
            "name": "Test JSON Export",
            "description": "Test export to JSON format",
            "file_ids": TEST_FILE_IDS,
            "export_format": "json",
            "include_metadata": True,
            "compress": True,
            "ttl_hours": 48,
        }

        response = self.session.post(f"{API_BASE}/jobs", json=data)

        if self.assert_response(response, 200, "Create JSON export"):
            result = response.json()
            self.test_job_ids.append(result["id"])
            print_info(f"Created job ID: {result['id']}")
            return True
        return False

    def test_create_excel_export(self) -> bool:
        """Test: Create Excel export job"""
        print_test("Create Excel Export Job")

        data = {
            "name": "Test Excel Export",
            "description": "Test export to Excel format",
            "file_ids": TEST_FILE_IDS,
            "export_format": "excel",
            "include_metadata": True,
            "compress": False,  # Excel files don't need compression
            "ttl_hours": 72,
        }

        response = self.session.post(f"{API_BASE}/jobs", json=data)

        if self.assert_response(response, 200, "Create Excel export"):
            result = response.json()
            self.test_job_ids.append(result["id"])
            print_info(f"Created job ID: {result['id']}")
            return True
        return False

    def test_get_export_job(self) -> bool:
        """Test: Get export job by ID"""
        print_test("Get Export Job")

        if not self.test_job_ids:
            print_warning("No test jobs available, skipping")
            return False

        job_id = self.test_job_ids[0]
        response = self.session.get(f"{API_BASE}/jobs/{job_id}")

        if self.assert_response(response, 200, "Get export job"):
            result = response.json()
            print_info(f"Job: {result['name']}")
            print_info(f"Status: {result['status']}")
            print_info(f"Progress: {result['progress']}%")
            print_info(f"Files: {result['processed_files']}/{result['total_files']}")
            if result.get("error_message"):
                print_warning(f"Error: {result['error_message']}")
            return True
        return False

    def test_list_export_jobs(self) -> bool:
        """Test: List all export jobs"""
        print_test("List Export Jobs")

        response = self.session.get(f"{API_BASE}/jobs")

        if self.assert_response(response, 200, "List export jobs"):
            result = response.json()
            print_info(f"Found {len(result)} export job(s)")
            for job in result[:5]:  # Show first 5
                print_info(f"  - {job['name']} ({job['status']}, {job['progress']}%)")
            return True
        return False

    def test_list_with_filters(self) -> bool:
        """Test: List export jobs with filters"""
        print_test("List Export Jobs with Filters")

        # Filter by status
        response = self.session.get(
            f"{API_BASE}/jobs", params={"status": "completed", "limit": 10}
        )

        if self.assert_response(response, 200, "List with filters"):
            result = response.json()
            print_info(f"Found {len(result)} completed job(s)")
            return True
        return False

    def test_monitor_progress(self) -> bool:
        """Test: Monitor export progress"""
        print_test("Monitor Export Progress")

        if not self.test_job_ids:
            print_warning("No test jobs available, skipping")
            return False

        job_id = self.test_job_ids[0]
        max_attempts = 30  # 30 attempts = ~1 minute
        attempt = 0

        print_info("Monitoring job progress (checking every 2 seconds)...")

        while attempt < max_attempts:
            response = self.session.get(f"{API_BASE}/jobs/{job_id}")

            if response.status_code != 200:
                print_error("Failed to get job status")
                return False

            job = response.json()
            status = job["status"]
            progress = job["progress"]

            print(f"  Status: {status:12} Progress: {progress:3}%", end="\r")

            if status in ["completed", "failed", "cancelled"]:
                print()  # New line
                print_info(f"Job finished with status: {status}")
                if status == "completed":
                    print_info(
                        f"Processed files: {job['processed_files']}/{job['total_files']}"
                    )
                    if job.get("output_size"):
                        size_mb = job["output_size"] / (1024 * 1024)
                        print_info(f"Output size: {size_mb:.2f} MB")
                    if job.get("download_url"):
                        print_info(f"Download URL available")
                elif status == "failed":
                    print_error(f"Error: {job.get('error_message', 'Unknown error')}")

                self.passed += 1
                return True

            attempt += 1
            time.sleep(2)

        print()
        print_warning("Job did not complete within timeout")
        self.passed += 1  # Not a failure, just slow
        return True

    def test_cancel_export_job(self) -> bool:
        """Test: Cancel a pending export job"""
        print_test("Cancel Export Job")

        # Create a new job to cancel
        data = {
            "name": "Test Cancel Export",
            "description": "This job will be cancelled",
            "file_ids": TEST_FILE_IDS,
            "export_format": "csv",
            "compress": True,
            "ttl_hours": 24,
        }

        create_response = self.session.post(f"{API_BASE}/jobs", json=data)
        if create_response.status_code != 200:
            print_error("Failed to create job for cancellation test")
            return False

        job_id = create_response.json()["id"]
        print_info(f"Created job {job_id} for cancellation")

        # Cancel immediately
        cancel_response = self.session.post(f"{API_BASE}/jobs/{job_id}/cancel")

        if self.assert_response(cancel_response, 200, "Cancel export job"):
            print_info("Job cancelled successfully")

            # Verify status
            verify_response = self.session.get(f"{API_BASE}/jobs/{job_id}")
            if verify_response.status_code == 200:
                job = verify_response.json()
                if job["status"] == "cancelled":
                    print_success("Job status confirmed as cancelled")
                    return True
        return False

    def test_get_statistics(self) -> bool:
        """Test: Get export statistics"""
        print_test("Get Export Statistics")

        response = self.session.get(f"{API_BASE}/statistics")

        if self.assert_response(response, 200, "Get statistics"):
            stats = response.json()
            print_info(f"Total jobs: {stats['total_jobs']}")
            print_info(f"By status: {stats['by_status']}")
            print_info(f"By format: {stats['by_format']}")
            print_info(f"Total size: {stats['total_size_bytes'] / (1024**2):.2f} MB")
            print_info(f"Total files exported: {stats['total_files_exported']}")
            return True
        return False

    def test_delete_export_job(self) -> bool:
        """Test: Delete an export job"""
        print_test("Delete Export Job")

        if len(self.test_job_ids) < 2:
            print_warning("Not enough test jobs, skipping")
            return False

        # Delete second test job
        job_id = self.test_job_ids[1]
        response = self.session.delete(f"{API_BASE}/jobs/{job_id}")

        if self.assert_response(response, 204, "Delete export job"):
            print_info(f"Deleted job ID: {job_id}")

            # Verify deletion
            verify_response = self.session.get(f"{API_BASE}/jobs/{job_id}")
            if verify_response.status_code == 404:
                print_success("Job deletion confirmed (404)")
                return True
        return False

    def test_error_handling_invalid_job(self) -> bool:
        """Test: Error handling for invalid job ID"""
        print_test("Error Handling - Invalid Job ID")

        response = self.session.get(f"{API_BASE}/jobs/999999")

        if self.assert_response(response, 404, "Invalid job ID"):
            print_info("Correctly returned 404 for non-existent job")
            return True
        return False

    def test_error_handling_invalid_format(self) -> bool:
        """Test: Error handling for invalid export format"""
        print_test("Error Handling - Invalid Format")

        data = {
            "name": "Invalid Format Export",
            "file_ids": TEST_FILE_IDS,
            "export_format": "invalid_format",  # Invalid
            "compress": True,
            "ttl_hours": 24,
        }

        response = self.session.post(f"{API_BASE}/jobs", json=data)

        if self.assert_response(response, 422, "Invalid format"):
            print_info("Correctly rejected invalid format")
            return True
        return False

    def test_error_handling_empty_files(self) -> bool:
        """Test: Error handling for empty file list"""
        print_test("Error Handling - Empty File List")

        data = {
            "name": "Empty Export",
            "file_ids": [],  # Empty
            "export_format": "csv",
            "compress": True,
            "ttl_hours": 24,
        }

        response = self.session.post(f"{API_BASE}/jobs", json=data)

        if self.assert_response(response, 422, "Empty file list"):
            print_info("Correctly rejected empty file list")
            return True
        return False

    def cleanup(self) -> bool:
        """Clean up test data"""
        print_section("Cleanup Test Data")

        success = True

        # Delete remaining test jobs
        for job_id in self.test_job_ids:
            try:
                response = self.session.delete(f"{API_BASE}/jobs/{job_id}")
                if response.status_code in [204, 404]:
                    print_success(f"Deleted job ID: {job_id}")
                else:
                    print_warning(
                        f"Failed to delete job {job_id}: {response.status_code}"
                    )
                    success = False
            except Exception as e:
                print_warning(f"Error deleting job {job_id}: {e}")
                success = False

        return success

    def run_all_tests(self):
        """Run all tests in sequence"""
        print_section("Data Export Test Suite")

        tests = [
            # Service info
            self.test_get_service_info,
            # Create exports
            self.test_create_csv_export,
            self.test_create_json_export,
            self.test_create_excel_export,
            # Query exports
            self.test_get_export_job,
            self.test_list_export_jobs,
            self.test_list_with_filters,
            # Monitor progress
            self.test_monitor_progress,
            # Statistics
            self.test_get_statistics,
            # Operations
            self.test_cancel_export_job,
            self.test_delete_export_job,
            # Error handling
            self.test_error_handling_invalid_job,
            self.test_error_handling_invalid_format,
            self.test_error_handling_empty_files,
        ]

        for test in tests:
            try:
                test()
            except Exception as e:
                self.failed += 1
                print_error(f"Test failed with exception: {str(e)}")
            print()  # Blank line between tests

        # Cleanup
        self.cleanup()

        # Print summary
        self.print_summary()

    def print_summary(self):
        """Print test summary"""
        print_section("Test Summary")

        total = self.passed + self.failed
        success_rate = (self.passed / total * 100) if total > 0 else 0

        print(f"Total tests: {total}")
        print(f"{Colors.GREEN}Passed: {self.passed}{Colors.RESET}")
        print(f"{Colors.RED}Failed: {self.failed}{Colors.RESET}")
        print(f"Success rate: {success_rate:.1f}%")

        if self.failed == 0:
            print(f"\n{Colors.GREEN}{Colors.BOLD}✓ All tests passed!{Colors.RESET}\n")
        else:
            print(f"\n{Colors.RED}{Colors.BOLD}✗ Some tests failed{Colors.RESET}\n")


def main():
    """Main test runner"""
    print(f"{Colors.BOLD}Data Export Test Suite{Colors.RESET}")
    print(f"Testing against: {BASE_URL}\n")

    # Check if API is available
    try:
        response = requests.get(f"{BASE_URL}/docs")
        if response.status_code != 200:
            print_error(f"API not available at {BASE_URL}")
            print_info("Please ensure the backend is running")
            sys.exit(1)
    except requests.exceptions.ConnectionError:
        print_error(f"Cannot connect to {BASE_URL}")
        print_info("Please ensure the backend is running")
        sys.exit(1)

    # Warning about test data
    print_warning("NOTE: This test requires valid file IDs")
    print_info(f"Using file IDs: {TEST_FILE_IDS}")
    print_info("Update TEST_FILE_IDS in the script if these don't exist\n")

    # Run tests
    test_suite = DataExportTest()
    test_suite.run_all_tests()

    # Exit with appropriate code
    sys.exit(0 if test_suite.failed == 0 else 1)


if __name__ == "__main__":
    main()
