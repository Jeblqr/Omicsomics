#!/usr/bin/env python3
"""
Test script to verify API response formats
Checks for consistency across all endpoints
"""

import requests
import json
from typing import Dict, Any

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
    }
    color = colors.get(status, Colors.RESET)
    print(f"{color}{msg}{Colors.RESET}")


def print_header(text: str):
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*70}")
    print(f"{text:^70}")
    print(f"{'='*70}{Colors.RESET}\n")


def print_json(data: Any):
    print(json.dumps(data, indent=2, default=str))


def login() -> str:
    """Login and get token"""
    response = requests.post(
        f"{BASE_URL}/auth/login/access-token",
        data={"username": "admin@example.com", "password": "admin123"},
    )
    if response.status_code == 200:
        return response.json()["access_token"]
    raise Exception("Login failed")


def test_projects_api(token: str):
    """Test projects CRUD operations"""
    print_header("Testing Projects API")
    headers = {"Authorization": f"Bearer {token}"}

    # 1. Create project
    print_status("1. Creating project...", "info")
    create_response = requests.post(
        f"{BASE_URL}/projects/",
        headers=headers,
        json={"name": "Test Project", "description": "Test description"},
    )
    print(f"Status: {create_response.status_code}")
    print(f"Response:")
    print_json(create_response.json())

    if create_response.status_code not in [200, 201]:
        print_status("❌ Create failed!", "error")
        return

    project = create_response.json()
    project_id = project.get("id")
    print_status(f"✅ Created project ID: {project_id}", "success")

    # 2. List projects
    print_status("\n2. Listing projects...", "info")
    list_response = requests.get(f"{BASE_URL}/projects/", headers=headers)
    print(f"Status: {list_response.status_code}")
    print(f"Response (first item):")
    if list_response.json():
        print_json(list_response.json()[0])
    else:
        print_status("⚠️  No projects found!", "warning")

    # 3. Get single project
    print_status(f"\n3. Getting project {project_id}...", "info")
    get_response = requests.get(f"{BASE_URL}/projects/{project_id}", headers=headers)
    print(f"Status: {get_response.status_code}")
    print(f"Response:")
    print_json(get_response.json())

    # 4. Update project
    print_status(f"\n4. Updating project {project_id}...", "info")
    update_response = requests.put(
        f"{BASE_URL}/projects/{project_id}",
        headers=headers,
        json={"name": "Updated Test Project", "description": "Updated description"},
    )
    print(f"Status: {update_response.status_code}")
    print(f"Response:")
    print_json(update_response.json())

    # 5. Delete project
    print_status(f"\n5. Deleting project {project_id}...", "info")
    delete_response = requests.delete(
        f"{BASE_URL}/projects/{project_id}", headers=headers
    )
    print(f"Status: {delete_response.status_code}")
    print(f"Response:")
    print_json(delete_response.json())
    print_status(
        f"✅ Delete response format: {type(delete_response.json())}", "success"
    )

    # 6. Verify deletion
    print_status(f"\n6. Verifying project {project_id} is deleted...", "info")
    verify_response = requests.get(f"{BASE_URL}/projects/{project_id}", headers=headers)
    print(f"Status: {verify_response.status_code}")
    if verify_response.status_code == 404:
        print_status("✅ Project successfully deleted!", "success")
    else:
        print_status("❌ Project still exists!", "error")


def test_runs_api(token: str, project_id: int):
    """Test runs CRUD operations"""
    print_header("Testing Runs API")
    headers = {"Authorization": f"Bearer {token}"}

    # 1. Create run
    print_status("1. Creating run...", "info")
    create_response = requests.post(
        f"{BASE_URL}/runs/",
        headers=headers,
        json={
            "name": "Test Run",
            "description": "Test description",
            "project_id": project_id,
            "pipeline_type": "template",
            "pipeline_template_id": "rna-seq-basic",
            "input_files": [],
            "parameters": {},
        },
    )
    print(f"Status: {create_response.status_code}")
    print(f"Response:")
    print_json(create_response.json())

    if create_response.status_code not in [200, 201]:
        print_status("❌ Create failed!", "error")
        return

    run = create_response.json()
    run_id = run.get("id")
    print_status(f"✅ Created run ID: {run_id}", "success")

    # 2. List runs
    print_status("\n2. Listing runs...", "info")
    list_response = requests.get(f"{BASE_URL}/runs/", headers=headers)
    print(f"Status: {list_response.status_code}")
    print(f"Response (first item):")
    if list_response.json():
        print_json(list_response.json()[0])

    # 3. Delete run
    print_status(f"\n3. Deleting run {run_id}...", "info")
    delete_response = requests.delete(f"{BASE_URL}/runs/{run_id}", headers=headers)
    print(f"Status: {delete_response.status_code}")
    print(f"Response:")
    print_json(delete_response.json())


def test_data_api(token: str, project_id: int):
    """Test data files CRUD operations"""
    print_header("Testing Data API")
    headers = {"Authorization": f"Bearer {token}"}

    # Create a test file
    import io

    test_file_content = b"sample,value\na,1\nb,2\n"
    files = {"file": ("test.csv", io.BytesIO(test_file_content), "text/csv")}
    data = {"project_id": str(project_id)}

    print_status("1. Uploading file...", "info")
    upload_response = requests.post(
        f"{BASE_URL}/data/upload", headers=headers, files=files, data=data
    )
    print(f"Status: {upload_response.status_code}")
    print(f"Response:")
    print_json(upload_response.json())

    if upload_response.status_code not in [200, 201]:
        print_status("❌ Upload failed!", "error")
        return

    file_data = upload_response.json()
    file_id = file_data.get("id")
    print_status(f"✅ Uploaded file ID: {file_id}", "success")

    # 2. List files
    print_status("\n2. Listing data files...", "info")
    list_response = requests.get(
        f"{BASE_URL}/data/?project_id={project_id}", headers=headers
    )
    print(f"Status: {list_response.status_code}")
    print(f"Response (first item):")
    if list_response.json():
        print_json(list_response.json()[0])

    # 3. Delete file
    print_status(f"\n3. Deleting file {file_id}...", "info")
    delete_response = requests.delete(f"{BASE_URL}/data/{file_id}", headers=headers)
    print(f"Status: {delete_response.status_code}")
    print(f"Response:")
    print_json(delete_response.json())


def main():
    print_header("🔍 API FORMAT VERIFICATION TEST 🔍")

    # Login
    print_status("Logging in...", "info")
    try:
        token = login()
        print_status("✅ Login successful\n", "success")
    except Exception as e:
        print_status(f"❌ Login failed: {e}", "error")
        return

    # Test projects
    test_projects_api(token)

    # Create a project for other tests
    headers = {"Authorization": f"Bearer {token}"}
    proj_response = requests.post(
        f"{BASE_URL}/projects/",
        headers=headers,
        json={"name": "Test Project for APIs", "description": "For testing"},
    )
    if proj_response.status_code in [200, 201]:
        project_id = proj_response.json()["id"]

        # Test runs
        test_runs_api(token, project_id)

        # Test data
        test_data_api(token, project_id)

        # Cleanup
        requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)

    print_header("TEST COMPLETE")
    print_status("\nCheck above for any inconsistencies in response formats", "info")


if __name__ == "__main__":
    main()
