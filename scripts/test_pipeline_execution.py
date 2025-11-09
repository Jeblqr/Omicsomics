#!/usr/bin/env python3
"""
Test script for pipeline execution functionality
Tests the complete flow: create run -> start execution -> check progress -> view logs
"""

import requests
import time
import json
from typing import Optional

BASE_URL = "http://localhost:8000"

class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

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

def login() -> Optional[str]:
    """Login and get access token"""
    print_header("Step 1: Login")
    
    try:
        response = requests.post(
            f"{BASE_URL}/auth/login/access-token",
            data={
                "username": "admin@example.com",
                "password": "admin123"
            }
        )
        
        if response.status_code == 200:
            token = response.json()["access_token"]
            print_status("✅ Login successful", "success")
            return token
        else:
            print_status(f"❌ Login failed: {response.status_code} - {response.text}", "error")
            return None
    except Exception as e:
        print_status(f"❌ Login error: {str(e)}", "error")
        return None

def get_projects(token: str) -> Optional[int]:
    """Get first project ID"""
    print_header("Step 2: Get Project")
    
    try:
        response = requests.get(
            f"{BASE_URL}/projects/",
            headers={"Authorization": f"Bearer {token}"}
        )
        
        if response.status_code == 200:
            projects = response.json()
            if projects:
                project_id = projects[0]["id"]
                project_name = projects[0]["name"]
                print_status(f"✅ Found project: {project_name} (ID: {project_id})", "success")
                return project_id
            else:
                print_status("⚠️  No projects found, creating one...", "warning")
                # Create a test project
                create_response = requests.post(
                    f"{BASE_URL}/projects/",
                    headers={"Authorization": f"Bearer {token}"},
                    json={
                        "name": "Pipeline Test Project",
                        "description": "Auto-created for pipeline testing"
                    }
                )
                if create_response.status_code == 201:
                    project_id = create_response.json()["id"]
                    print_status(f"✅ Created test project (ID: {project_id})", "success")
                    return project_id
                else:
                    print_status(f"❌ Failed to create project: {create_response.text}", "error")
                    return None
        else:
            print_status(f"❌ Failed to get projects: {response.status_code}", "error")
            return None
    except Exception as e:
        print_status(f"❌ Error: {str(e)}", "error")
        return None

def create_run(token: str, project_id: int) -> Optional[int]:
    """Create a new pipeline run"""
    print_header("Step 3: Create Pipeline Run")
    
    try:
        run_data = {
            "name": "Test RNA-seq Pipeline Run",
            "description": "Testing pipeline executor with RNA-seq template",
            "project_id": project_id,
            "pipeline_type": "template",
            "pipeline_template_id": "rna-seq-basic",
            "input_files": [],
            "parameters": {
                "threads": 8,
                "min_quality": 30
            }
        }
        
        print_status(f"Creating run with template: {run_data['pipeline_template_id']}", "info")
        
        response = requests.post(
            f"{BASE_URL}/runs/",
            headers={"Authorization": f"Bearer {token}"},
            json=run_data
        )
        
        if response.status_code == 201:
            run = response.json()
            run_id = run["id"]
            print_status(f"✅ Run created successfully (ID: {run_id})", "success")
            print_status(f"   Name: {run['name']}", "info")
            print_status(f"   Status: {run['status']}", "info")
            print_status(f"   Pipeline: {run['pipeline_template_id']}", "info")
            return run_id
        else:
            print_status(f"❌ Failed to create run: {response.status_code} - {response.text}", "error")
            return None
    except Exception as e:
        print_status(f"❌ Error: {str(e)}", "error")
        return None

def start_run(token: str, run_id: int) -> bool:
    """Start pipeline execution"""
    print_header("Step 4: Start Pipeline Execution")
    
    try:
        response = requests.post(
            f"{BASE_URL}/runs/{run_id}/start",
            headers={"Authorization": f"Bearer {token}"}
        )
        
        if response.status_code == 200:
            print_status("✅ Pipeline execution started!", "success")
            print_status("   Run is now executing in the background...", "info")
            return True
        else:
            print_status(f"❌ Failed to start run: {response.status_code} - {response.text}", "error")
            return False
    except Exception as e:
        print_status(f"❌ Error: {str(e)}", "error")
        return False

def monitor_run(token: str, run_id: int, max_wait: int = 120):
    """Monitor run progress until completion"""
    print_header("Step 5: Monitor Execution Progress")
    
    start_time = time.time()
    last_progress = -1
    last_status = ""
    
    print_status("Polling run status every 2 seconds...\n", "info")
    
    try:
        while time.time() - start_time < max_wait:
            response = requests.get(
                f"{BASE_URL}/runs/{run_id}",
                headers={"Authorization": f"Bearer {token}"}
            )
            
            if response.status_code == 200:
                run = response.json()
                status = run["status"]
                progress = run.get("progress", 0)
                
                # Print update if changed
                if progress != last_progress or status != last_status:
                    elapsed = int(time.time() - start_time)
                    
                    if status == "running":
                        bar_length = 40
                        filled = int(bar_length * progress / 100)
                        bar = "█" * filled + "░" * (bar_length - filled)
                        print_status(
                            f"[{elapsed:3d}s] Status: {status:12s} | Progress: {bar} {progress:5.1f}%",
                            "cyan"
                        )
                    else:
                        status_color = "success" if status == "completed" else "error" if status == "failed" else "warning"
                        print_status(
                            f"[{elapsed:3d}s] Status: {status:12s} | Progress: {progress:5.1f}%",
                            status_color
                        )
                    
                    last_progress = progress
                    last_status = status
                
                # Check if finished
                if status in ["completed", "failed", "cancelled"]:
                    print()
                    if status == "completed":
                        print_status("✅ Pipeline execution COMPLETED!", "success")
                    elif status == "failed":
                        print_status("❌ Pipeline execution FAILED!", "error")
                        if run.get("error_message"):
                            print_status(f"   Error: {run['error_message']}", "error")
                    else:
                        print_status("⚠️  Pipeline execution CANCELLED", "warning")
                    
                    return run
            else:
                print_status(f"⚠️  Failed to get run status: {response.status_code}", "warning")
            
            time.sleep(2)
        
        print_status(f"\n⏱️  Timeout after {max_wait}s", "warning")
        return None
        
    except Exception as e:
        print_status(f"\n❌ Monitoring error: {str(e)}", "error")
        return None

def view_logs(token: str, run_id: int):
    """View run logs"""
    print_header("Step 6: View Execution Logs")
    
    try:
        response = requests.get(
            f"{BASE_URL}/runs/{run_id}/logs",
            headers={"Authorization": f"Bearer {token}"}
        )
        
        if response.status_code == 200:
            data = response.json()
            logs = data.get("logs", "")
            error_msg = data.get("error_message", "")
            
            if logs:
                print_status("Logs:", "cyan")
                print(f"{Colors.RESET}{logs}")
            else:
                print_status("No logs available", "warning")
            
            if error_msg:
                print()
                print_status("Error Message:", "error")
                print(f"{Colors.RED}{error_msg}{Colors.RESET}")
            
            return True
        else:
            print_status(f"❌ Failed to get logs: {response.status_code}", "error")
            return False
    except Exception as e:
        print_status(f"❌ Error: {str(e)}", "error")
        return False

def main():
    print_header("🧬 PIPELINE EXECUTION TEST 🧬")
    
    # Step 1: Login
    token = login()
    if not token:
        print_status("\n❌ Test FAILED: Cannot login", "error")
        return
    
    # Step 2: Get project
    project_id = get_projects(token)
    if not project_id:
        print_status("\n❌ Test FAILED: Cannot get/create project", "error")
        return
    
    # Step 3: Create run
    run_id = create_run(token, project_id)
    if not run_id:
        print_status("\n❌ Test FAILED: Cannot create run", "error")
        return
    
    # Step 4: Start execution
    if not start_run(token, run_id):
        print_status("\n❌ Test FAILED: Cannot start run", "error")
        return
    
    # Step 5: Monitor progress
    final_run = monitor_run(token, run_id, max_wait=180)  # 3 minutes max
    
    # Step 6: View logs
    view_logs(token, run_id)
    
    # Final result
    print_header("TEST RESULT")
    if final_run and final_run["status"] == "completed":
        print_status("🎉 ✅ ALL TESTS PASSED! Pipeline execution works perfectly!", "success")
        print_status(f"   Final Progress: {final_run.get('progress', 0)}%", "success")
        print_status(f"   Total Steps: {len([l for l in final_run.get('logs', '').split('\\n') if 'Step' in l])}", "success")
    elif final_run and final_run["status"] == "failed":
        print_status("⚠️  Test completed but pipeline FAILED", "warning")
        print_status("   Check logs above for error details", "warning")
    else:
        print_status("⚠️  Test incomplete or timeout", "warning")
    
    print()

if __name__ == "__main__":
    main()
