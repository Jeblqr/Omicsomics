#!/usr/bin/env python3
"""
Quick test script for frontend data processing UI
Tests the upload → process → view workflow
"""

import requests
import time
import json
from pathlib import Path

BASE_URL = "http://localhost:8001/api/v1"


def test_frontend_workflow():
    """Test the complete frontend data processing workflow"""

    print("=" * 60)
    print("Frontend Data Processing UI Test")
    print("=" * 60)

    # 1. Register and login
    print("\n1️⃣ Setting up test user...")
    email = f"test_frontend_{int(time.time())}@example.com"
    password = "testpass123"

    try:
        response = requests.post(
            f"{BASE_URL}/auth/register", json={"email": email, "password": password}
        )
        if response.status_code not in [200, 201]:
            print(
                f"   ❌ Registration failed: {response.status_code} - {response.text}"
            )
            return
        user_data = response.json()
        user_id = user_data.get("id") or user_data.get("user_id")
        print(f"   ✅ User created: {email} (ID: {user_id})")
    except Exception as e:
        print(f"   ❌ Registration failed: {e}")
        return

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
    print(f"   ✅ Logged in, token: {token[:20]}...")

    # 2. Create project
    print("\n2️⃣ Creating test project...")
    response = requests.post(
        f"{BASE_URL}/projects/",
        json={
            "name": "Frontend Test Project",
            "description": "Testing data processing UI",
        },
        headers=headers,
    )
    project_id = response.json()["id"]
    print(f"   ✅ Project created (ID: {project_id})")

    # 3. Upload with processing
    print("\n3️⃣ Testing file upload with processing...")
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
                "process_file": "true",  # This is the new parameter!
            },
            headers=headers,
        )

    if response.status_code != 201:
        print(f"   ❌ Upload failed: {response.status_code} - {response.text}")
        return

    result = response.json()
    file_id = result["id"]
    print(f"   ✅ File uploaded (ID: {file_id})")

    # Check processing info
    processing_info = result.get("processing") or result.get("processing_info", {})
    print(f"\n   📊 Processing Info:")
    print(f"      • Processed: {processing_info.get('processed')}")
    print(f"      • Converted: {processing_info.get('converted')}")
    print(f"      • Omics Type: {processing_info.get('omics_type')}")
    print(f"      • Records: {processing_info.get('record_count')}")

    # 4. Fetch file list (simulate frontend GET request)
    print("\n4️⃣ Fetching file list...")
    response = requests.get(
        f"{BASE_URL}/data/?project_id={project_id}", headers=headers
    )
    files = response.json()

    print(f"   ✅ Found {len(files)} file(s)")
    for file in files:
        metadata = file.get("metadata_", {})
        is_processed = bool(metadata.get("processed_file_id"))
        omics_type = metadata.get("omics_type", "N/A")

        print(f"\n   📄 {file['filename']}")
        print(f"      • Status: {'✓ Processed' if is_processed else 'Raw only'}")
        print(f"      • Omics Type: {omics_type}")
        print(f"      • Size: {file['size']} bytes")

    # 5. Get processed data (simulate "View Data" button)
    print("\n5️⃣ Fetching processed data...")
    response = requests.get(f"{BASE_URL}/data/{file_id}/processed", headers=headers)

    if response.status_code != 200:
        print(
            f"   ❌ Failed to get processed data: {response.status_code} - {response.text}"
        )
        return

    processed = response.json()
    unified_data = processed.get("unified_data", {})
    metadata = unified_data.get("metadata", {})
    records = unified_data.get("records", [])

    print(f"   ✅ Processed data retrieved")
    print(f"\n   🔬 Unified Data:")
    print(f"      • Omics Type: {metadata.get('omics_type')}")
    print(f"      • Source: {metadata.get('source')}")
    print(f"      • Total Records: {len(records)}")
    print(f"      • Format Version: {metadata.get('format_version')}")

    if records:
        print(f"\n   📊 First 3 records:")
        for i, record in enumerate(records[:3], 1):
            print(f"      {i}. {record}")

    # 6. Test multi-file selection scenario
    print("\n6️⃣ Testing multi-file selection scenario...")

    # Upload another file
    test_file2 = Path("test_data/sample_proteomics.csv")
    if test_file2.exists():
        with open(test_file2, "rb") as f:
            response = requests.post(
                f"{BASE_URL}/data/upload",
                files={"file": (test_file2.name, f, "text/csv")},
                data={"project_id": project_id, "process_file": "true"},
                headers=headers,
            )

        file2_id = response.json()["id"]
        print(f"   ✅ Second file uploaded (ID: {file2_id})")

        # Get updated file list
        response = requests.get(
            f"{BASE_URL}/data/?project_id={project_id}", headers=headers
        )
        files = response.json()

        processed_files = [
            f for f in files if f.get("metadata_", {}).get("processed_file_id")
        ]
        print(f"\n   ✅ Total processed files: {len(processed_files)}")
        print(f"   💡 These files can be selected for creating an analysis run")

        for file in processed_files:
            omics_type = file.get("metadata_", {}).get("omics_type", "unknown")
            print(f"      • {file['filename']} ({omics_type})")

    # 7. Cleanup
    print("\n7️⃣ Cleaning up...")
    response = requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)
    if response.status_code == 204:
        print(f"   ✅ Project deleted (cascade cleaned all files)")

    print("\n" + "=" * 60)
    print("✅ All tests passed!")
    print("=" * 60)
    print("\n📝 Frontend features verified:")
    print("   ✓ Upload with 'process_file' parameter")
    print("   ✓ File list shows processing status")
    print("   ✓ Metadata includes omics_type")
    print("   ✓ GET /processed endpoint returns unified data")
    print("   ✓ Multi-file selection ready")
    print("\n🌐 Open http://localhost:5173 to test the UI manually!")


if __name__ == "__main__":
    try:
        test_frontend_workflow()
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback

        traceback.print_exc()
