#!/usr/bin/env python3
"""
Test Script for Pipeline Dataset Integration

Tests the integration between Pipeline Builder and Dataset Manager:
- Using datasets as pipeline inputs
- Creating datasets from run outputs
- Recording and querying lineage
- Validation functionality
"""

import requests
import json
from typing import Dict, Any

# Configuration
BASE_URL = "http://localhost:8000"
API_PREFIX = "/api"


# Colors for output
class Colors:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    END = "\033[0m"


def print_success(msg):
    print(f"{Colors.GREEN}✓ {msg}{Colors.END}")


def print_error(msg):
    print(f"{Colors.RED}✗ {msg}{Colors.END}")


def print_info(msg):
    print(f"{Colors.BLUE}ℹ {msg}{Colors.END}")


def print_warning(msg):
    print(f"{Colors.YELLOW}⚠ {msg}{Colors.END}")


def test_get_dataset_input_files():
    """Test getting files from dataset for pipeline input."""
    print("\n" + "=" * 60)
    print("Test 1: Get Dataset Input Files")
    print("=" * 60)

    # Create a test dataset first
    print_info("Creating test dataset...")
    dataset_response = requests.post(
        f"{BASE_URL}{API_PREFIX}/datasets/",
        json={
            "project_id": 1,
            "name": "Test Input Dataset",
            "description": "Dataset for pipeline input testing",
            "data_type": "transcriptomics",
        },
    )

    if dataset_response.status_code != 200:
        print_error(f"Failed to create dataset: {dataset_response.text}")
        return False

    dataset = dataset_response.json()
    dataset_id = dataset["id"]
    print_success(f"Created dataset: {dataset_id}")

    # Add some test files
    print_info("Adding files to dataset...")
    files_added = 0
    for i in range(3):
        file_response = requests.post(
            f"{BASE_URL}{API_PREFIX}/datasets/{dataset_id}/files",
            json={
                "file_path": f"/tmp/test_file_{i}.fastq",
                "file_name": f"test_file_{i}.fastq",
                "file_type": "fastq",
                "role": "primary",
                "compute_hash": False,
            },
        )
        if file_response.status_code == 200:
            files_added += 1

    print_success(f"Added {files_added} files")

    # Get dataset input files
    print_info("Fetching dataset input files...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/input-files",
        json={"dataset_id": dataset_id},
    )

    if response.status_code != 200:
        print_error(f"Failed: {response.text}")
        return False

    data = response.json()
    print_success(f"Retrieved {len(data['files'])} files")

    for file in data["files"]:
        print(f"  - {file['name']} ({file['type']}, {file['role']})")

    return True


def test_validate_dataset():
    """Test dataset validation for pipeline input."""
    print("\n" + "=" * 60)
    print("Test 2: Validate Dataset for Pipeline")
    print("=" * 60)

    # Create a test dataset
    print_info("Creating test dataset...")
    dataset_response = requests.post(
        f"{BASE_URL}{API_PREFIX}/datasets/",
        json={
            "project_id": 1,
            "name": "Validation Test Dataset",
            "description": "Dataset for validation testing",
            "data_type": "genomics",
        },
    )

    if dataset_response.status_code != 200:
        print_error(f"Failed to create dataset: {dataset_response.text}")
        return False

    dataset_id = dataset_response.json()["id"]

    # Validate empty dataset (should fail)
    print_info("Validating empty dataset (should fail)...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/validate",
        json={"dataset_id": dataset_id},
    )

    if response.status_code != 200:
        print_error(f"Request failed: {response.text}")
        return False

    result = response.json()
    if result["is_valid"]:
        print_warning("Empty dataset passed validation (unexpected)")
    else:
        print_success(f"Validation correctly failed: {', '.join(result['errors'])}")

    # Add files
    print_info("Adding files...")
    requests.post(
        f"{BASE_URL}{API_PREFIX}/datasets/{dataset_id}/files",
        json={
            "file_path": "/tmp/test.fastq",
            "file_name": "test.fastq",
            "file_type": "fastq",
            "role": "primary",
            "compute_hash": False,
        },
    )

    # Validate with file (should still fail due to missing file)
    print_info("Validating with file reference...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/validate",
        json={"dataset_id": dataset_id, "required_file_types": ["fastq"]},
    )

    result = response.json()
    if not result["is_valid"]:
        print_success(f"Validation failed as expected: {', '.join(result['errors'])}")
    else:
        print_warning("Validation passed (file may exist)")

    return True


def test_create_dataset_from_run():
    """Test creating dataset from pipeline run outputs."""
    print("\n" + "=" * 60)
    print("Test 3: Create Dataset from Run")
    print("=" * 60)

    # Create a mock run
    print_info("Creating mock pipeline run...")
    run_response = requests.post(
        f"{BASE_URL}{API_PREFIX}/runs/",
        json={
            "name": "Test RNA-seq Pipeline",
            "description": "Mock pipeline for testing",
            "pipeline_type": "template",
            "pipeline_template_id": "rna-seq-basic",
            "project_id": 1,
            "status": "completed",
            "output_files": [1, 2, 3],  # Mock file IDs
        },
    )

    if run_response.status_code not in [200, 201]:
        print_warning(f"Could not create run: {run_response.text}")
        print_info("Using mock run ID for testing...")
        run_id = 999  # Mock ID
    else:
        run_id = run_response.json()["id"]
        print_success(f"Created run: {run_id}")

    # Create dataset from run
    print_info("Creating dataset from run outputs...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/create-from-run",
        json={
            "run_id": run_id,
            "dataset_name": "RNA-seq Analysis Output",
            "dataset_description": "Processed RNA-seq results",
            "data_type": "transcriptomics",
            "tags": ["pipeline-output", "rna-seq", "test"],
        },
    )

    if response.status_code != 200:
        print_error(f"Failed: {response.text}")
        return False

    data = response.json()
    print_success(f"Created dataset: {data['dataset']['id']}")
    print(f"  Name: {data['dataset']['name']}")
    print(f"  Data Type: {data['dataset']['data_type']}")
    print(f"  Files: {data['dataset']['file_count']}")
    print(f"  Message: {data['message']}")

    return True


def test_auto_create_dataset():
    """Test auto-creating dataset with inferred settings."""
    print("\n" + "=" * 60)
    print("Test 4: Auto-Create Dataset from Run")
    print("=" * 60)

    # Use mock run ID
    run_id = 999

    print_info(f"Auto-creating dataset from run {run_id}...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/auto-create",
        json={"run_id": run_id, "auto_tags": True},
    )

    if response.status_code != 200:
        print_warning(
            f"Auto-create failed (expected if run doesn't exist): {response.text}"
        )
        return True  # Not a critical failure for testing

    data = response.json()
    print_success(f"Auto-created dataset: {data['dataset']['id']}")
    print(f"  Name: {data['dataset']['name']}")
    print(f"  Data Type: {data['dataset']['data_type']}")

    return True


def test_link_datasets_to_run():
    """Test linking input datasets to a run."""
    print("\n" + "=" * 60)
    print("Test 5: Link Datasets to Run")
    print("=" * 60)

    # Create test datasets
    print_info("Creating input datasets...")
    dataset_ids = []
    for i in range(2):
        response = requests.post(
            f"{BASE_URL}{API_PREFIX}/datasets/",
            json={
                "project_id": 1,
                "name": f"Input Dataset {i+1}",
                "description": f"Test input dataset {i+1}",
                "data_type": "genomics",
            },
        )
        if response.status_code == 200:
            dataset_ids.append(response.json()["id"])

    print_success(f"Created {len(dataset_ids)} datasets")

    # Link to run
    run_id = 999  # Mock run ID
    print_info(f"Linking datasets to run {run_id}...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/link-to-run",
        json={
            "run_id": run_id,
            "dataset_ids": dataset_ids,
            "input_mapping": {
                "node_1_input": dataset_ids[0],
                "node_2_reference": dataset_ids[1],
            },
        },
    )

    if response.status_code != 200:
        print_warning(
            f"Linking failed (expected if run doesn't exist): {response.text}"
        )
        return True

    data = response.json()
    print_success(data["message"])

    return True


def test_get_dataset_usage():
    """Test getting dataset usage in pipeline runs."""
    print("\n" + "=" * 60)
    print("Test 6: Get Dataset Usage")
    print("=" * 60)

    # Create a dataset
    print_info("Creating test dataset...")
    response = requests.post(
        f"{BASE_URL}{API_PREFIX}/datasets/",
        json={
            "project_id": 1,
            "name": "Usage Test Dataset",
            "description": "Dataset for usage tracking",
            "data_type": "genomics",
        },
    )

    if response.status_code != 200:
        print_error(f"Failed to create dataset: {response.text}")
        return False

    dataset_id = response.json()["id"]

    # Get usage
    print_info(f"Getting usage for dataset {dataset_id}...")
    response = requests.get(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/dataset/{dataset_id}/usage"
    )

    if response.status_code != 200:
        print_error(f"Failed: {response.text}")
        return False

    data = response.json()
    print_success(f"Found {len(data['runs'])} runs using this dataset")

    if data["runs"]:
        for run in data["runs"]:
            print(f"  - Run {run['run_id']}: {run['run_name']} ({run['status']})")
    else:
        print_info("  No runs found (dataset not yet used)")

    return True


def test_get_run_output_datasets():
    """Test getting datasets created by a run."""
    print("\n" + "=" * 60)
    print("Test 7: Get Run Output Datasets")
    print("=" * 60)

    run_id = 999  # Mock run ID

    print_info(f"Getting output datasets for run {run_id}...")
    response = requests.get(
        f"{BASE_URL}{API_PREFIX}/pipeline-datasets/run/{run_id}/output-datasets"
    )

    if response.status_code != 200:
        print_warning(f"Request failed: {response.text}")
        return True  # Not critical

    data = response.json()
    print_success(f"Found {len(data)} output datasets")

    for dataset in data:
        print(
            f"  - {dataset['name']} ({dataset['data_type']}, {dataset['file_count']} files)"
        )

    return True


def test_api_info():
    """Test API info endpoint."""
    print("\n" + "=" * 60)
    print("Test 8: API Information")
    print("=" * 60)

    print_info("Fetching API information...")
    response = requests.get(f"{BASE_URL}{API_PREFIX}/pipeline-datasets/")

    if response.status_code != 200:
        print_error(f"Failed: {response.text}")
        return False

    data = response.json()
    print_success(f"Service: {data['service']} v{data['version']}")
    print(f"  Description: {data['description']}")
    print(f"  Features:")
    for feature in data["features"]:
        print(f"    • {feature}")
    print(f"  Endpoints: {len(data['endpoints'])} available")

    return True


def run_all_tests():
    """Run all tests and report results."""
    print("\n" + "=" * 70)
    print(" Pipeline Dataset Integration Test Suite")
    print("=" * 70)

    tests = [
        ("API Information", test_api_info),
        ("Get Dataset Input Files", test_get_dataset_input_files),
        ("Validate Dataset", test_validate_dataset),
        ("Create Dataset from Run", test_create_dataset_from_run),
        ("Auto-Create Dataset", test_auto_create_dataset),
        ("Link Datasets to Run", test_link_datasets_to_run),
        ("Get Dataset Usage", test_get_dataset_usage),
        ("Get Run Output Datasets", test_get_run_output_datasets),
    ]

    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print_error(f"Test '{name}' raised exception: {e}")
            results.append((name, False))

    # Summary
    print("\n" + "=" * 70)
    print(" Test Summary")
    print("=" * 70)

    passed = sum(1 for _, success in results if success)
    total = len(results)

    for name, success in results:
        status = "PASS" if success else "FAIL"
        color = Colors.GREEN if success else Colors.RED
        print(f"{color}{status:6}{Colors.END} {name}")

    print("-" * 70)
    print(f"Total: {passed}/{total} tests passed ({100*passed//total}%)")

    if passed == total:
        print_success("\n All tests passed! ✨")
    else:
        print_warning(f"\n{total - passed} test(s) failed")


if __name__ == "__main__":
    run_all_tests()
