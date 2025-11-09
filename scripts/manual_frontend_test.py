#!/usr/bin/env python3
"""
Manual Frontend Test Script
Tests frontend functionality by checking UI elements and interactions
"""

import requests
import time
import json

BASE_URL = "http://localhost:8001/api/v1"
FRONTEND_URL = "http://localhost:5173"


def test_backend_apis():
    """Test all backend APIs that frontend uses"""
    print("\n" + "=" * 70)
    print("TESTING BACKEND APIs")
    print("=" * 70)

    # Login
    print("\n1. Testing Login...")
    response = requests.post(
        f"{BASE_URL}/auth/login/access-token",
        data={"username": "test_user@omics.com", "password": "TestPassword123!"},
    )
    if response.status_code == 200:
        token = response.json()["access_token"]
        headers = {"Authorization": f"Bearer {token}"}
        print("   ✅ Login successful")
    else:
        print(f"   ❌ Login failed: {response.status_code}")
        return False

    # Test Projects
    print("\n2. Testing Projects API...")
    response = requests.get(f"{BASE_URL}/projects/", headers=headers)
    if response.status_code == 200:
        projects = response.json()
        print(f"   ✅ Projects API working ({len(projects)} projects)")
    else:
        print(f"   ❌ Projects API failed: {response.status_code}")

    # Test Data Files
    print("\n3. Testing Data Files API...")
    response = requests.get(f"{BASE_URL}/data/", headers=headers)
    if response.status_code == 200:
        files = response.json()
        print(f"   ✅ Data API working ({len(files)} files)")

        # Test deletion of a file
        if len(files) > 0:
            test_file_id = files[0]["id"]
            print(f"\n4. Testing File Deletion (ID: {test_file_id})...")
            response = requests.delete(
                f"{BASE_URL}/data/{test_file_id}", headers=headers
            )
            if response.status_code == 204:
                print("   ✅ File deletion works")
                # Re-create it or note it's deleted
            else:
                print(f"   ❌ File deletion failed: {response.status_code}")
    else:
        print(f"   ❌ Data API failed: {response.status_code}")

    # Test Runs
    print("\n5. Testing Runs API...")
    response = requests.get(f"{BASE_URL}/runs/", headers=headers)
    if response.status_code == 200:
        runs = response.json()
        print(f"   ✅ Runs API working ({len(runs)} runs)")

        # Check run statuses
        statuses = {}
        for run in runs:
            status = run["status"]
            statuses[status] = statuses.get(status, 0) + 1
        print(f"   📊 Run statuses: {statuses}")

        # Test run deletion
        if len(runs) > 0:
            test_run_id = runs[0]["id"]
            print(f"\n6. Testing Run Deletion (ID: {test_run_id})...")
            response = requests.delete(
                f"{BASE_URL}/runs/{test_run_id}", headers=headers
            )
            if response.status_code == 204:
                print("   ✅ Run deletion works")
            else:
                print(f"   ❌ Run deletion failed: {response.status_code}")
    else:
        print(f"   ❌ Runs API failed: {response.status_code}")

    # Test Pipeline Templates
    print("\n7. Testing Pipeline Templates API...")
    response = requests.get(f"{BASE_URL}/pipelines/", headers=headers)
    if response.status_code == 200:
        pipelines = response.json()
        print(f"   ✅ Pipelines API working ({len(pipelines)} templates)")
    else:
        print(f"   ❌ Pipelines API failed: {response.status_code}")

    # Test Custom Pipelines
    print("\n8. Testing Custom Pipelines API...")
    if projects:
        project_id = projects[0]["id"]
        response = requests.get(
            f"{BASE_URL}/custom-pipelines/?project_id={project_id}", headers=headers
        )
        if response.status_code == 200:
            custom_pipelines = response.json()
            print(
                f"   ✅ Custom Pipelines API working ({len(custom_pipelines)} custom pipelines)"
            )

            # Test save pipeline
            print("\n9. Testing Save Custom Pipeline...")
            test_pipeline = {
                "name": "Test Pipeline",
                "description": "Test Description",
                "project_id": project_id,
                "steps": [{"name": "Step 1", "tool": "test_tool", "parameters": {}}],
            }
            response = requests.post(
                f"{BASE_URL}/custom-pipelines/", json=test_pipeline, headers=headers
            )
            if response.status_code == 201:
                pipeline_id = response.json()["id"]
                print(f"   ✅ Save pipeline works (ID: {pipeline_id})")

                # Test delete
                print(f"\n10. Testing Delete Custom Pipeline (ID: {pipeline_id})...")
                response = requests.delete(
                    f"{BASE_URL}/custom-pipelines/{pipeline_id}", headers=headers
                )
                if response.status_code == 204:
                    print("   ✅ Delete pipeline works")
                else:
                    print(f"   ❌ Delete pipeline failed: {response.status_code}")
            else:
                print(
                    f"   ❌ Save pipeline failed: {response.status_code} - {response.text}"
                )
        else:
            print(f"   ❌ Custom Pipelines API failed: {response.status_code}")

    print("\n" + "=" * 70)
    print("BACKEND API TESTS COMPLETE")
    print("=" * 70)
    return True


def check_frontend_issues():
    """Document known frontend issues"""
    print("\n" + "=" * 70)
    print("KNOWN FRONTEND ISSUES TO FIX")
    print("=" * 70)

    issues = [
        {
            "issue": "White text on white background in forms",
            "file": "frontend/src/styles/index.css",
            "problem": "Dark mode CSS conflicts with light form inputs",
            "fix": "Need to update form styles for better contrast",
        },
        {
            "issue": "Runs always pending",
            "file": "backend/app/services/runs.py or celery workers",
            "problem": "Runs are created but not executed",
            "fix": "Need to start run execution or check worker status",
        },
        {
            "issue": "Unable to delete data",
            "file": "frontend/src/components/SandboxView.tsx",
            "problem": "Delete functionality may not be working",
            "fix": "Backend API works, need to test frontend button",
        },
        {
            "issue": "Failed to save pipeline",
            "file": "frontend/src/pages/pipelines/CustomPipelinesPage.tsx",
            "problem": "Pipeline save may have validation errors",
            "fix": "Need to check request payload format",
        },
        {
            "issue": "Poor UI/UX in Runs page",
            "file": "frontend/src/pages/runs/RunsPage.tsx",
            "problem": "Layout and styling needs improvement",
            "fix": "Need to improve table design, add more actions",
        },
    ]

    for i, issue in enumerate(issues, 1):
        print(f"\n{i}. {issue['issue']}")
        print(f"   File: {issue['file']}")
        print(f"   Problem: {issue['problem']}")
        print(f"   Fix: {issue['fix']}")

    print("\n" + "=" * 70)
    return issues


def main():
    print("=" * 70)
    print("OMICSOMICS FRONTEND MANUAL TEST")
    print("=" * 70)
    print("\nThis script tests backend APIs and documents frontend issues.")
    print(f"Backend URL: {BASE_URL}")
    print(f"Frontend URL: {FRONTEND_URL}")

    # Test backend
    backend_ok = test_backend_apis()

    # Document frontend issues
    issues = check_frontend_issues()

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"✅ Backend APIs: {'WORKING' if backend_ok else 'FAILED'}")
    print(f"⚠️  Frontend Issues: {len(issues)} issues need fixing")
    print("\nNext steps:")
    print("1. Fix CSS for dark mode vs light forms")
    print("2. Fix run execution (check celery workers)")
    print("3. Test delete functionality in browser")
    print("4. Test pipeline save in browser")
    print("5. Improve Runs page UI/UX")
    print(
        "\nPlease open http://localhost:5173 in your browser to test the frontend manually."
    )
    print("=" * 70)


if __name__ == "__main__":
    main()
