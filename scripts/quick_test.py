#!/usr/bin/env python3
"""
Quick API Health Check and Basic Workflow Test
"""

import requests
import json
import sys
from pathlib import Path

BASE_URL = "http://localhost:8001/api/v1"


def test_health():
    """Test health endpoint"""
    print("Testing health endpoint...")
    try:
        response = requests.get(f"{BASE_URL}/health/")
        print(f"✅ Health check: {response.json()}")
        return True
    except Exception as e:
        print(f"❌ Health check failed: {e}")
        return False


def test_auth():
    """Test authentication"""
    print("\nTesting authentication...")
    try:
        # Register
        reg_data = {
            "email": "quicktest@example.com",
            "password": "testpass123",
            "full_name": "Quick Test User",
        }
        requests.post(f"{BASE_URL}/auth/register", json=reg_data)

        # Login
        login_data = {"username": "quicktest@example.com", "password": "testpass123"}
        response = requests.post(f"{BASE_URL}/auth/login", data=login_data)
        token = response.json()["access_token"]
        print(f"✅ Authentication successful")
        return token
    except Exception as e:
        print(f"❌ Authentication failed: {e}")
        return None


def test_project_crud(token):
    """Test project CRUD operations"""
    print("\nTesting project operations...")
    headers = {"Authorization": f"Bearer {token}"}

    try:
        # Create
        project = {"name": "Quick Test Project", "description": "Test"}
        response = requests.post(f"{BASE_URL}/projects/", json=project, headers=headers)
        project_id = response.json()["id"]
        print(f"✅ Created project ID: {project_id}")

        # Read
        response = requests.get(f"{BASE_URL}/projects/", headers=headers)
        projects = response.json()
        print(f"✅ Retrieved {len(projects)} projects")

        # Delete
        requests.delete(f"{BASE_URL}/projects/{project_id}", headers=headers)
        print(f"✅ Deleted project ID: {project_id}")

        return True
    except Exception as e:
        print(f"❌ Project operations failed: {e}")
        return False


def test_pipeline_list(token):
    """Test pipeline listing"""
    print("\nTesting pipeline templates...")
    headers = {"Authorization": f"Bearer {token}"}

    try:
        response = requests.get(f"{BASE_URL}/pipelines/", headers=headers)
        pipelines = response.json()
        print(f"✅ Found {len(pipelines)} pipeline templates:")
        for p in pipelines[:5]:
            print(f"   - {p['id']}: {p['name']}")
        return True
    except Exception as e:
        print(f"❌ Pipeline listing failed: {e}")
        return False


def main():
    print("=" * 60)
    print("QUICK OMICS PLATFORM TEST")
    print("=" * 60)

    results = []

    # Run tests
    results.append(("Health Check", test_health()))

    token = test_auth()
    if token:
        results.append(("Authentication", True))
        results.append(("Project CRUD", test_project_crud(token)))
        results.append(("Pipeline List", test_pipeline_list(token)))
    else:
        results.append(("Authentication", False))

    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    for name, result in results:
        status = "✅ PASSED" if result else "❌ FAILED"
        print(f"{name}: {status}")

    passed = sum(1 for _, r in results if r)
    total = len(results)
    print(f"\n{passed}/{total} tests passed ({passed/total*100:.0f}%)")

    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
