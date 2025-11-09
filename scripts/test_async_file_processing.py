#!/usr/bin/env python3
"""
Test script for async file processing with Celery
Tests the full async flow: upload → submit task → poll status → retrieve result
"""

import asyncio
import httpx
import time
import sys
from pathlib import Path

# Test configuration
BASE_URL = "http://localhost:8001"
TEST_FILE = Path(__file__).parent.parent / "test_data" / "sample_transcriptomics.csv"


async def main():
    print("=== Async File Processing Test ===\n")

    # Create a test CSV file if it doesn't exist
    if not TEST_FILE.exists():
        TEST_FILE.parent.mkdir(parents=True, exist_ok=True)
        TEST_FILE.write_text(
            "gene_id,sample1,sample2,sample3\n"
            "ENSG00000000003,100,200,150\n"
            "ENSG00000000005,50,75,60\n"
            "ENSG00000000419,300,400,350\n"
        )
        print(f"✓ Created test file: {TEST_FILE}\n")

    async with httpx.AsyncClient(timeout=30.0) as client:
        # Step 1: Health check
        print("1. Checking backend health...")
        try:
            response = await client.get(f"{BASE_URL}/healthz")
            if response.status_code == 200:
                print(f"✓ Backend is healthy: {response.json()}\n")
            else:
                print(f"✗ Backend health check failed: {response.status_code}")
                return
        except Exception as e:
            print(f"✗ Cannot connect to backend: {e}")
            print("Make sure backend and celery worker are running:\n")
            print("  cd backend && uvicorn app.main:app --reload")
            print("  celery -A app.celery_app worker --loglevel=info\n")
            return

        # Step 2: Create a test project (or use existing)
        print("2. Creating test project...")
        project_data = {
            "name": "Async Processing Test",
            "description": "Test project for async file processing",
        }

        # Note: This assumes authentication is handled
        # In real scenario, you'd need to login first
        try:
            response = await client.post(f"{BASE_URL}/projects", json=project_data)
            if response.status_code in [200, 201]:
                project = response.json()
                project_id = project["id"]
                print(f"✓ Created project: {project['name']} (ID: {project_id})\n")
            else:
                # Try to use existing project
                response = await client.get(f"{BASE_URL}/projects")
                if response.status_code == 200:
                    projects = response.json()
                    if projects:
                        project_id = projects[0]["id"]
                        print(f"✓ Using existing project (ID: {project_id})\n")
                    else:
                        print("✗ No projects available. Please create one manually.")
                        return
                else:
                    print(f"✗ Failed to create/list projects: {response.status_code}")
                    return
        except Exception as e:
            print(f"✗ Project creation failed: {e}")
            return

        # Step 3: Upload file with async processing
        print("3. Uploading file with async processing enabled...")
        with open(TEST_FILE, "rb") as f:
            files = {"file": (TEST_FILE.name, f, "text/csv")}
            data = {
                "project_id": str(project_id),
                "process_file": "true",
                "async_processing": "true",  # Enable async processing
                "sample_id": "test_sample_001",
            }

            try:
                response = await client.post(
                    f"{BASE_URL}/data/upload", files=files, data=data
                )

                if response.status_code in [200, 201]:
                    result = response.json()
                    print(f"✓ File uploaded successfully!")
                    print(f"  File ID: {result['id']}")
                    print(f"  Filename: {result['filename']}")

                    processing_info = result.get("processing", {})
                    if processing_info.get("async"):
                        task_id = processing_info.get("task_id")
                        print(f"  Task ID: {task_id}")
                        print(f"  Status: {processing_info.get('status')}\n")

                        # Step 4: Poll task status
                        print("4. Polling task status...")
                        max_wait = 60  # seconds
                        poll_interval = 2  # seconds
                        elapsed = 0

                        while elapsed < max_wait:
                            status_response = await client.get(
                                f"{BASE_URL}/data/task/{task_id}/status"
                            )

                            if status_response.status_code == 200:
                                status = status_response.json()
                                state = status.get("state")
                                progress = status.get("progress", 0)
                                message = status.get("status", "Processing...")

                                print(
                                    f"  [{elapsed}s] State: {state}, Progress: {progress}% - {message}"
                                )

                                if state == "SUCCESS":
                                    print(f"\n✓ Processing completed successfully!")
                                    result_data = status.get("result", {})
                                    print(
                                        f"  Processed file ID: {result_data.get('processed_file_id')}"
                                    )
                                    print(
                                        f"  Omics type: {result_data.get('omics_type')}"
                                    )
                                    print(
                                        f"  Record count: {result_data.get('record_count')}"
                                    )
                                    break
                                elif state == "FAILURE":
                                    print(f"\n✗ Processing failed!")
                                    print(f"  Error: {status.get('error')}")
                                    break
                                elif state in ["PENDING", "STARTED", "PROGRESS"]:
                                    # Continue polling
                                    await asyncio.sleep(poll_interval)
                                    elapsed += poll_interval
                                else:
                                    print(f"\n✗ Unknown state: {state}")
                                    break
                            else:
                                print(
                                    f"✗ Failed to get task status: {status_response.status_code}"
                                )
                                break

                        if elapsed >= max_wait:
                            print(
                                f"\n⚠ Timeout: Task did not complete in {max_wait} seconds"
                            )

                    else:
                        print("⚠ Async processing was not enabled in response")
                        if processing_info.get("error"):
                            print(f"  Error: {processing_info.get('error')}")

                else:
                    print(f"✗ Upload failed: {response.status_code}")
                    print(f"  Response: {response.text}")
                    return

            except Exception as e:
                print(f"✗ Upload failed: {e}")
                return

        print("\n=== Test Complete ===")


if __name__ == "__main__":
    asyncio.run(main())
