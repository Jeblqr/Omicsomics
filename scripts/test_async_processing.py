#!/usr/bin/env python3
"""
Test script for async file processing with Celery
Tests the upload → async process → check status workflow
"""

import requests
import time
import json
from pathlib import Path

BASE_URL = "http://localhost:8001/api/v1"


def test_async_processing():
    """Test async file processing workflow"""

    print("=" * 60)
    print("Async File Processing Test")
    print("=" * 60)

    # 1. Setup user
    print("\n1️⃣ Setting up test user...")
    email = f"test_async_{int(time.time())}@example.com"
    password = "testpass123"

    response = requests.post(
        f"{BASE_URL}/auth/register", json={"email": email, "password": password}
    )

    if response.status_code not in [200, 201]:
        print(f"   ❌ Registration failed: {response.status_code} - {response.text}")
        return

    user_id = response.json().get("id")
    print(f"   ✅ User created (ID: {user_id})")

    # Login
    response = requests.post(
        f"{BASE_URL}/auth/login/access-token",
        data={"username": email, "password": password},
    )

    if response.status_code != 200:
        print(f"   ❌ Login failed: {response.status_code} - {response.text}")
        return

    token = response.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}
    print(f"   ✅ Logged in")

    # 2. Create project
    print("\n2️⃣ Creating project...")
    response = requests.post(
        f"{BASE_URL}/projects/",
        json={
            "name": "Async Processing Test",
            "description": "Testing Celery async file processing",
        },
        headers=headers,
    )

    project_id = response.json()["id"]
    print(f"   ✅ Project created (ID: {project_id})")

    # 3. Upload file with async processing
    print("\n3️⃣ Uploading file with async processing...")
    test_file = Path("test_data/sample_transcriptomics.csv")

    if not test_file.exists():
        print(f"   ❌ Test file not found: {test_file}")
        return

    with open(test_file, "rb") as f:
        response = requests.post(
            f"{BASE_URL}/data/upload",
            files={"file": (test_file.name, f, "text/csv")},
            data={
                "project_id": project_id,
                "process_file": "true",
                "async_processing": "true",  # KEY: Enable async processing
            },
            headers=headers,
        )

    if response.status_code != 201:
        print(f"   ❌ Upload failed: {response.status_code} - {response.text}")
        return

    result = response.json()
    file_id = result["id"]
    processing = result.get("processing", {})

    print(f"   ✅ File uploaded (ID: {file_id})")
    print(f"      • Async: {processing.get('async')}")
    print(f"      • Task ID: {processing.get('task_id')}")
    print(f"      • Status: {processing.get('status')}")

    task_id = processing.get("task_id")
    if not task_id:
        print("   ❌ No task_id returned")
        return

    # 4. Poll task status
    print("\n4️⃣ Monitoring task progress...")
    max_polls = 30
    poll_interval = 2

    for i in range(max_polls):
        response = requests.get(
            f"{BASE_URL}/data/task/{task_id}/status", headers=headers
        )

        if response.status_code != 200:
            print(f"   ❌ Status check failed: {response.status_code}")
            break

        status_data = response.json()
        state = status_data.get("state")
        status_msg = status_data.get("status")
        progress = status_data.get("progress", 0)

        print(
            f"   [{i+1:2d}/{max_polls}] State: {state:10s} Progress: {progress:3d}% - {status_msg}"
        )

        if state == "SUCCESS":
            result = status_data.get("result", {})
            print(f"\n   ✅ Processing complete!")
            print(f"      • Processed File ID: {result.get('processed_file_id')}")
            print(f"      • Omics Type: {result.get('omics_type')}")
            print(f"      • Record Count: {result.get('record_count')}")

            # 5. Verify processed data
            print("\n5️⃣ Verifying processed data...")
            response = requests.get(
                f"{BASE_URL}/data/{file_id}/processed", headers=headers
            )

            if response.status_code == 200:
                processed = response.json()
                unified_data = processed.get("unified_data", {})
                records = unified_data.get("records", [])
                print(f"   ✅ Retrieved {len(records)} records")
                print(f"      • First record: {records[0] if records else 'N/A'}")
            else:
                print(f"   ⚠️ Could not retrieve processed data: {response.status_code}")

            break

        elif state == "FAILURE":
            error = status_data.get("error")
            print(f"\n   ❌ Processing failed: {error}")
            break

        time.sleep(poll_interval)

    else:
        print(f"\n   ⚠️ Timeout waiting for task completion")

    # 6. Cleanup
    print("\n6️⃣ Cleaning up...")
    response = requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)
    if response.status_code == 204:
        print(f"   ✅ Project deleted")

    print("\n" + "=" * 60)
    print("Test Complete!")
    print("=" * 60)
    print("\n📝 Key features tested:")
    print("   ✓ Upload with async_processing=true")
    print("   ✓ Receive task_id immediately")
    print("   ✓ Poll task status endpoint")
    print("   ✓ Monitor progress updates")
    print("   ✓ Retrieve processed data after completion")


if __name__ == "__main__":
    try:
        test_async_processing()
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback

        traceback.print_exc()
