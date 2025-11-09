#!/usr/bin/env python3
"""
Test script to verify:
1. Projects can be deleted (with cascade delete of runs/data/samples)
2. Runs progress updates correctly (not stuck at 0%)
"""

import requests
import time
import json

BASE_URL = "http://localhost:8000"


class Colors:
    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    RESET = "\033[0m"
    BOLD = "\033[1m"


def print_status(msg: str, status: str = "info"):
    colors = {
        "success": Colors.GREEN,
        "error": Colors.RED,
        "warning": Colors.YELLOW,
        "info": Colors.BLUE,
        "cyan": Colors.CYAN,
    }
    color = colors.get(status, Colors.RESET)
    print(f"{color}{msg}{Colors.RESET}")


def print_header(text: str):
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*70}")
    print(f"{text:^70}")
    print(f"{'='*70}{Colors.RESET}\n")


def login():
    response = requests.post(
        f"{BASE_URL}/auth/login/access-token",
        data={"username": "admin@example.com", "password": "admin123"},
    )
    if response.status_code == 200:
        return response.json()["access_token"]
    raise Exception("Login failed")


def test_project_deletion_with_dependencies(token: str):
    """Test that projects can be deleted even when they have runs and data files"""
    print_header("Test 1: Project Deletion with Dependencies")
    headers = {"Authorization": f"Bearer {token}"}

    # Step 1: Create a project
    print_status("Step 1: Creating project...", "info")
    proj_resp = requests.post(
        f"{BASE_URL}/projects/",
        headers=headers,
        json={"name": "Delete Test Project", "description": "Will be deleted"},
    )
    if proj_resp.status_code != 201:
        print_status(f"❌ Failed to create project: {proj_resp.status_code}", "error")
        return False
    project_id = proj_resp.json()["id"]
    print_status(f"✅ Created project ID: {project_id}", "success")

    # Step 2: Upload a data file to this project
    print_status("\nStep 2: Uploading data file...", "info")
    import io

    files = {"file": ("test.csv", io.BytesIO(b"sample,value\na,1\nb,2\n"), "text/csv")}
    data = {"project_id": str(project_id)}
    upload_resp = requests.post(
        f"{BASE_URL}/data/upload", headers=headers, files=files, data=data
    )
    if upload_resp.status_code != 201:
        print_status(
            f"⚠️  Warning: Failed to upload file: {upload_resp.status_code}", "warning"
        )
    else:
        print_status(f"✅ Uploaded data file", "success")

    # Step 3: Create a run in this project
    print_status("\nStep 3: Creating run...", "info")
    run_resp = requests.post(
        f"{BASE_URL}/runs/",
        headers=headers,
        json={
            "name": "Test Run",
            "description": "Test",
            "project_id": project_id,
            "pipeline_type": "template",
            "pipeline_template_id": "rna-seq-basic",
            "input_files": [],
            "parameters": {},
        },
    )
    if run_resp.status_code != 201:
        print_status(
            f"⚠️  Warning: Failed to create run: {run_resp.status_code}", "warning"
        )
    else:
        print_status(f"✅ Created run", "success")

    # Step 4: Try to delete the project
    print_status("\nStep 4: Deleting project (with all dependencies)...", "info")
    delete_resp = requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)

    if delete_resp.status_code == 204:
        print_status("✅ Project deleted successfully!", "success")

        # Step 5: Verify project is gone
        print_status("\nStep 5: Verifying project is deleted...", "info")
        verify_resp = requests.get(f"{BASE_URL}/projects/{project_id}", headers=headers)
        if verify_resp.status_code == 404:
            print_status("✅ Confirmed: Project no longer exists", "success")
            return True
        else:
            print_status("❌ Error: Project still exists!", "error")
            return False
    else:
        print_status(f"❌ Delete failed with status {delete_resp.status_code}", "error")
        print_status(f"Response: {delete_resp.text}", "error")
        return False


def test_run_progress_updates(token: str):
    """Test that run progress updates correctly (not stuck at 0%)"""
    print_header("Test 2: Run Progress Updates")
    headers = {"Authorization": f"Bearer {token}"}

    # Step 1: Create a project
    print_status("Step 1: Creating project...", "info")
    proj_resp = requests.post(
        f"{BASE_URL}/projects/",
        headers=headers,
        json={"name": "Progress Test", "description": "Test"},
    )
    if proj_resp.status_code != 201:
        print_status(f"❌ Failed to create project", "error")
        return False
    project_id = proj_resp.json()["id"]
    print_status(f"✅ Created project ID: {project_id}", "success")

    # Step 2: Create a run
    print_status("\nStep 2: Creating run...", "info")
    run_resp = requests.post(
        f"{BASE_URL}/runs/",
        headers=headers,
        json={
            "name": "Progress Test Run",
            "description": "Testing progress updates",
            "project_id": project_id,
            "pipeline_type": "template",
            "pipeline_template_id": "rna-seq-basic",
            "input_files": [],
            "parameters": {},
        },
    )
    if run_resp.status_code != 201:
        print_status(f"❌ Failed to create run", "error")
        requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)
        return False
    run_id = run_resp.json()["id"]
    print_status(f"✅ Created run ID: {run_id}", "success")

    # Step 3: Start the run
    print_status("\nStep 3: Starting run...", "info")
    start_resp = requests.post(f"{BASE_URL}/runs/{run_id}/start", headers=headers)
    if start_resp.status_code != 200:
        print_status(f"❌ Failed to start run", "error")
        requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)
        return False
    print_status("✅ Run started", "success")

    # Step 4: Monitor progress
    print_status("\nStep 4: Monitoring progress (20 seconds)...", "info")
    progress_values = []
    for i in range(10):
        time.sleep(2)
        check_resp = requests.get(f"{BASE_URL}/runs/{run_id}", headers=headers)
        if check_resp.status_code == 200:
            run_data = check_resp.json()
            progress = run_data.get("progress", 0)
            status = run_data.get("status", "unknown")
            progress_values.append(progress)

            bar_length = 30
            filled = int(bar_length * progress / 100)
            bar = "█" * filled + "░" * (bar_length - filled)
            print_status(
                f"  [{i*2:2d}s] Status: {status:12s} | Progress: {bar} {progress:5.1f}%",
                "cyan",
            )

            if status in ["completed", "failed", "cancelled"]:
                break

    # Step 5: Check if progress changed
    print_status("\nStep 5: Analyzing progress...", "info")
    print_status(f"Progress values: {progress_values}", "info")

    # Cleanup
    requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)

    # Determine result
    if len(set(progress_values)) > 1 and max(progress_values) > 0:
        print_status("✅ SUCCESS: Progress updated correctly!", "success")
        print_status(
            f"   Progress changed from {min(progress_values)}% to {max(progress_values)}%",
            "success",
        )
        return True
    else:
        print_status("❌ FAIL: Progress stuck at 0%", "error")
        return False


def main():
    print_header("🔧 FIX VERIFICATION TEST 🔧")

    # Login
    print_status("Logging in...", "info")
    try:
        token = login()
        print_status("✅ Login successful\n", "success")
    except Exception as e:
        print_status(f"❌ Login failed: {e}", "error")
        return

    # Run tests
    test1_passed = test_project_deletion_with_dependencies(token)
    test2_passed = test_run_progress_updates(token)

    # Final result
    print_header("FINAL RESULT")
    print_status(
        f"Test 1 (Project Deletion): {'✅ PASSED' if test1_passed else '❌ FAILED'}",
        "success" if test1_passed else "error",
    )
    print_status(
        f"Test 2 (Run Progress): {'✅ PASSED' if test2_passed else '❌ FAILED'}",
        "success" if test2_passed else "error",
    )

    if test1_passed and test2_passed:
        print_status("\n🎉 ALL TESTS PASSED! Both issues are fixed!", "success")
    else:
        print_status("\n⚠️  Some tests failed. Check logs above for details.", "warning")


if __name__ == "__main__":
    main()
