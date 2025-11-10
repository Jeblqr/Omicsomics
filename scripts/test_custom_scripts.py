#!/usr/bin/env python3
"""
Comprehensive Test Script for Custom Script Tools

Tests all functionality of the custom script management and execution system.
"""

import json
import os
import sys
import time
from typing import Dict, Any, Optional

import requests
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)

# Configuration
API_BASE_URL = os.getenv("API_URL", "http://localhost:8000")
API_URL = f"{API_BASE_URL}/api/custom-scripts"

# Test data
TEST_PYTHON_SCRIPT = """import pandas as pd
import json
import sys
import os

def main(params):
    # Simple test script
    input_value = params.get('input_value', 'default')
    multiplier = params.get('multiplier', 2)
    
    result = {
        'input': input_value,
        'output': input_value * multiplier,
        'status': 'success'
    }
    
    # Save output
    os.makedirs('output', exist_ok=True)
    with open('output/result.json', 'w') as f:
        json.dump(result, f)
    
    return result

if __name__ == '__main__':
    with open(sys.argv[2], 'r') as f:
        params = json.load(f)
    
    result = main(params)
    print(json.dumps(result))
"""

TEST_R_SCRIPT = """library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
params <- fromJSON(args[2])

# Simple calculation
result <- list(
  input = params$input_value,
  doubled = params$input_value * 2,
  status = "success"
)

# Save output
dir.create('output', showWarnings = FALSE)
write(toJSON(result), 'output/result.json')

cat(toJSON(result))
"""

TEST_BASH_SCRIPT = """#!/bin/bash

PARAMS_FILE=$2
INPUT=$(jq -r '.input_value // "default"' $PARAMS_FILE)

mkdir -p output
echo "Input: $INPUT" > output/result.txt

echo "{\\"status\\": \\"success\\", \\"input\\": \\"$INPUT\\"}"
"""


class CustomScriptTester:
    """Test harness for Custom Script Tools."""

    def __init__(self, token: Optional[str] = None):
        self.token = token or os.getenv("API_TOKEN")
        self.headers = {"Authorization": f"Bearer {self.token}"} if self.token else {}
        self.test_results = []
        self.script_id = None
        self.execution_id = None

    def print_header(self, text: str):
        """Print a formatted header."""
        print(f"\n{Fore.CYAN}{'=' * 80}")
        print(f"{Fore.CYAN}{text:^80}")
        print(f"{Fore.CYAN}{'=' * 80}{Style.RESET_ALL}\n")

    def print_success(self, text: str):
        """Print a success message."""
        print(f"{Fore.GREEN}✓ {text}{Style.RESET_ALL}")

    def print_error(self, text: str):
        """Print an error message."""
        print(f"{Fore.RED}✗ {text}{Style.RESET_ALL}")

    def print_info(self, text: str):
        """Print an info message."""
        print(f"{Fore.YELLOW}ℹ {text}{Style.RESET_ALL}")

    def print_result(self, data: Dict[str, Any]):
        """Print formatted result data."""
        print(f"{Fore.BLUE}{json.dumps(data, indent=2)}{Style.RESET_ALL}")

    def make_request(
        self, method: str, endpoint: str, **kwargs
    ) -> Optional[requests.Response]:
        """Make an API request."""
        url = f"{API_URL}{endpoint}"
        try:
            response = requests.request(
                method, url, headers=self.headers, timeout=30, **kwargs
            )
            return response
        except Exception as e:
            self.print_error(f"Request failed: {e}")
            return None

    def test_create_python_script(self) -> bool:
        """Test creating a Python script."""
        self.print_header("Test 1: Create Python Script")

        payload = {
            "name": "Test Python Script",
            "description": "Simple Python test script",
            "language": "python",
            "script_content": TEST_PYTHON_SCRIPT,
            "parameters_schema": {
                "type": "object",
                "properties": {
                    "input_value": {"type": "string", "default": "hello"},
                    "multiplier": {"type": "integer", "default": 2},
                },
            },
            "requirements": ["pandas"],
            "timeout": 300,
            "max_memory": 512,
            "visibility": "private",
            "category": "testing",
            "tags": ["test", "python"],
        }

        response = self.make_request("POST", "/scripts", json=payload)

        if response and response.status_code == 201:
            data = response.json()
            self.script_id = data["id"]
            self.print_success(f"Created script: {data['name']} (ID: {data['id']})")
            self.print_result(data)
            return True
        else:
            self.print_error(
                f"Failed to create script: {response.status_code if response else 'No response'}"
            )
            if response:
                print(response.text)
            return False

    def test_get_script(self) -> bool:
        """Test getting a script by ID."""
        self.print_header("Test 2: Get Script Details")

        if not self.script_id:
            self.print_error("No script ID available")
            return False

        response = self.make_request("GET", f"/scripts/{self.script_id}")

        if response and response.status_code == 200:
            data = response.json()
            self.print_success(f"Retrieved script: {data['name']}")
            self.print_info(f"Language: {data['language']}")
            self.print_info(f"Total executions: {data['total_executions']}")
            return True
        else:
            self.print_error(
                f"Failed to get script: {response.status_code if response else 'No response'}"
            )
            return False

    def test_list_scripts(self) -> bool:
        """Test listing scripts with filters."""
        self.print_header("Test 3: List Scripts")

        # Test without filters
        response = self.make_request("GET", "/scripts")

        if response and response.status_code == 200:
            data = response.json()
            self.print_success(
                f"Listed {len(data['items'])} scripts (total: {data['total']})"
            )

            # Test with language filter
            response = self.make_request("GET", "/scripts?language=python")
            if response and response.status_code == 200:
                data = response.json()
                self.print_success(f"Filtered by Python: {len(data['items'])} scripts")

            # Test with search
            response = self.make_request("GET", "/scripts?search=test")
            if response and response.status_code == 200:
                data = response.json()
                self.print_success(f"Searched 'test': {len(data['items'])} scripts")

            return True
        else:
            self.print_error(
                f"Failed to list scripts: {response.status_code if response else 'No response'}"
            )
            return False

    def test_validate_parameters(self) -> bool:
        """Test parameter validation."""
        self.print_header("Test 4: Validate Parameters")

        if not self.script_id:
            self.print_error("No script ID available")
            return False

        # Test valid parameters
        response = self.make_request(
            "POST",
            f"/scripts/{self.script_id}/validate-parameters",
            json={"parameters": {"input_value": "test", "multiplier": 3}},
        )

        if response and response.status_code == 200:
            data = response.json()
            if data["valid"]:
                self.print_success("Valid parameters accepted")
            else:
                self.print_error(f"Parameters rejected: {data['errors']}")
                return False

            # Test invalid parameters (wrong type)
            response = self.make_request(
                "POST",
                f"/scripts/{self.script_id}/validate-parameters",
                json={"parameters": {"input_value": 123, "multiplier": "invalid"}},
            )

            if response and response.status_code == 200:
                data = response.json()
                if not data["valid"]:
                    self.print_success("Invalid parameters correctly rejected")
                else:
                    self.print_error("Invalid parameters were accepted")
                    return False

            return True
        else:
            self.print_error(
                f"Failed to validate parameters: {response.status_code if response else 'No response'}"
            )
            return False

    def test_execute_script(self) -> bool:
        """Test executing a script."""
        self.print_header("Test 5: Execute Script")

        if not self.script_id:
            self.print_error("No script ID available")
            return False

        payload = {
            "parameters": {"input_value": "hello", "multiplier": 3},
            "description": "Test execution",
        }

        response = self.make_request(
            "POST", f"/scripts/{self.script_id}/execute", json=payload
        )

        if response and response.status_code == 202:
            data = response.json()
            self.execution_id = data["id"]
            self.print_success(
                f"Execution started: {data['execution_key']} (ID: {data['id']})"
            )
            self.print_info(f"Status: {data['status']}")
            return True
        else:
            self.print_error(
                f"Failed to execute script: {response.status_code if response else 'No response'}"
            )
            if response:
                print(response.text)
            return False

    def test_monitor_execution(self) -> bool:
        """Test monitoring execution progress."""
        self.print_header("Test 6: Monitor Execution Progress")

        if not self.execution_id:
            self.print_error("No execution ID available")
            return False

        max_attempts = 30
        attempt = 0

        while attempt < max_attempts:
            response = self.make_request("GET", f"/executions/{self.execution_id}")

            if not response or response.status_code != 200:
                self.print_error(
                    f"Failed to get execution: {response.status_code if response else 'No response'}"
                )
                return False

            data = response.json()
            status = data["status"]

            self.print_info(f"Attempt {attempt + 1}: Status = {status}")

            if status == "completed":
                self.print_success("Execution completed successfully")
                self.print_info(f"Duration: {data.get('duration', 'N/A')}s")
                self.print_info(f"Exit code: {data.get('exit_code', 'N/A')}")
                if data.get("output_text"):
                    self.print_info(f"Output: {data['output_text'][:200]}")
                return True
            elif status == "failed":
                self.print_error("Execution failed")
                if data.get("error_text"):
                    print(f"Error: {data['error_text']}")
                return False
            elif status in ["pending", "running"]:
                time.sleep(2)
                attempt += 1
            else:
                self.print_error(f"Unknown status: {status}")
                return False

        self.print_error("Execution timeout")
        return False

    def test_list_executions(self) -> bool:
        """Test listing executions."""
        self.print_header("Test 7: List Executions")

        response = self.make_request("GET", "/executions")

        if response and response.status_code == 200:
            data = response.json()
            self.print_success(
                f"Listed {len(data['items'])} executions (total: {data['total']})"
            )

            # Test with script filter
            if self.script_id:
                response = self.make_request(
                    "GET", f"/executions?script_id={self.script_id}"
                )
                if response and response.status_code == 200:
                    data = response.json()
                    self.print_success(
                        f"Filtered by script: {len(data['items'])} executions"
                    )

            return True
        else:
            self.print_error(
                f"Failed to list executions: {response.status_code if response else 'No response'}"
            )
            return False

    def test_get_templates(self) -> bool:
        """Test getting script templates."""
        self.print_header("Test 8: Get Script Templates")

        response = self.make_request("GET", "/templates")

        if response and response.status_code == 200:
            templates = response.json()
            self.print_success(f"Retrieved {len(templates)} templates")
            for template in templates:
                self.print_info(f"- {template['language']}: {template['name']}")
            return True
        else:
            self.print_error(
                f"Failed to get templates: {response.status_code if response else 'No response'}"
            )
            return False

    def test_get_languages(self) -> bool:
        """Test getting supported languages."""
        self.print_header("Test 9: Get Supported Languages")

        response = self.make_request("GET", "/languages")

        if response and response.status_code == 200:
            languages = response.json()
            self.print_success(f"Supported languages: {', '.join(languages)}")
            return True
        else:
            self.print_error(
                f"Failed to get languages: {response.status_code if response else 'No response'}"
            )
            return False

    def test_update_script(self) -> bool:
        """Test updating a script."""
        self.print_header("Test 10: Update Script")

        if not self.script_id:
            self.print_error("No script ID available")
            return False

        payload = {
            "description": "Updated test script description",
            "tags": ["test", "python", "updated"],
        }

        response = self.make_request("PUT", f"/scripts/{self.script_id}", json=payload)

        if response and response.status_code == 200:
            data = response.json()
            self.print_success("Script updated successfully")
            self.print_info(f"New description: {data['description']}")
            self.print_info(f"New tags: {', '.join(data['tags'])}")
            return True
        else:
            self.print_error(
                f"Failed to update script: {response.status_code if response else 'No response'}"
            )
            return False

    def test_create_r_script(self) -> bool:
        """Test creating an R script."""
        self.print_header("Test 11: Create R Script")

        payload = {
            "name": "Test R Script",
            "description": "Simple R test script",
            "language": "r",
            "script_content": TEST_R_SCRIPT,
            "parameters_schema": {
                "type": "object",
                "properties": {"input_value": {"type": "number", "default": 10}},
            },
            "requirements": ["jsonlite"],
            "timeout": 300,
            "visibility": "private",
            "category": "testing",
            "tags": ["test", "r"],
        }

        response = self.make_request("POST", "/scripts", json=payload)

        if response and response.status_code == 201:
            data = response.json()
            self.print_success(f"Created R script: {data['name']} (ID: {data['id']})")
            return True
        else:
            self.print_error(
                f"Failed to create R script: {response.status_code if response else 'No response'}"
            )
            return False

    def test_create_bash_script(self) -> bool:
        """Test creating a Bash script."""
        self.print_header("Test 12: Create Bash Script")

        payload = {
            "name": "Test Bash Script",
            "description": "Simple Bash test script",
            "language": "bash",
            "script_content": TEST_BASH_SCRIPT,
            "parameters_schema": {
                "type": "object",
                "properties": {"input_value": {"type": "string", "default": "test"}},
            },
            "requirements": ["jq"],
            "timeout": 300,
            "visibility": "private",
            "category": "testing",
            "tags": ["test", "bash"],
        }

        response = self.make_request("POST", "/scripts", json=payload)

        if response and response.status_code == 201:
            data = response.json()
            self.print_success(
                f"Created Bash script: {data['name']} (ID: {data['id']})"
            )
            return True
        else:
            self.print_error(
                f"Failed to create Bash script: {response.status_code if response else 'No response'}"
            )
            return False

    def test_delete_script(self) -> bool:
        """Test deleting a script."""
        self.print_header("Test 13: Delete Script")

        if not self.script_id:
            self.print_error("No script ID available")
            return False

        response = self.make_request("DELETE", f"/scripts/{self.script_id}")

        if response and response.status_code == 204:
            self.print_success(f"Script {self.script_id} deleted successfully")
            return True
        else:
            self.print_error(
                f"Failed to delete script: {response.status_code if response else 'No response'}"
            )
            return False

    def test_error_handling(self) -> bool:
        """Test error handling."""
        self.print_header("Test 14: Error Handling")

        # Test getting non-existent script
        response = self.make_request("GET", "/scripts/99999")
        if response and response.status_code == 404:
            self.print_success("404 error handled correctly for non-existent script")
        else:
            self.print_error("Failed to handle non-existent script")
            return False

        # Test invalid script creation
        response = self.make_request("POST", "/scripts", json={"name": "Test"})
        if response and response.status_code in [400, 422]:
            self.print_success("400/422 error handled correctly for invalid data")
        else:
            self.print_error("Failed to handle invalid script data")
            return False

        return True

    def run_all_tests(self):
        """Run all tests."""
        print(f"\n{Fore.MAGENTA}{'=' * 80}")
        print(f"{Fore.MAGENTA}Custom Script Tools - Comprehensive Test Suite")
        print(f"{Fore.MAGENTA}{'=' * 80}{Style.RESET_ALL}\n")

        self.print_info(f"API URL: {API_URL}")
        self.print_info(f"Token: {'Set' if self.token else 'Not set'}\n")

        tests = [
            ("Create Python Script", self.test_create_python_script),
            ("Get Script Details", self.test_get_script),
            ("List Scripts", self.test_list_scripts),
            ("Validate Parameters", self.test_validate_parameters),
            ("Execute Script", self.test_execute_script),
            ("Monitor Execution", self.test_monitor_execution),
            ("List Executions", self.test_list_executions),
            ("Get Templates", self.test_get_templates),
            ("Get Languages", self.test_get_languages),
            ("Update Script", self.test_update_script),
            ("Create R Script", self.test_create_r_script),
            ("Create Bash Script", self.test_create_bash_script),
            ("Delete Script", self.test_delete_script),
            ("Error Handling", self.test_error_handling),
        ]

        results = []
        for name, test_func in tests:
            try:
                result = test_func()
                results.append((name, result))
            except Exception as e:
                self.print_error(f"Test '{name}' raised exception: {e}")
                results.append((name, False))

        # Print summary
        self.print_header("Test Summary")
        passed = sum(1 for _, result in results if result)
        total = len(results)

        for name, result in results:
            if result:
                self.print_success(f"{name}")
            else:
                self.print_error(f"{name}")

        print(f"\n{Fore.CYAN}{'=' * 80}")
        print(f"{Fore.CYAN}Results: {passed}/{total} tests passed")
        print(f"{Fore.CYAN}{'=' * 80}{Style.RESET_ALL}\n")

        return passed == total


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Test Custom Script Tools functionality"
    )
    parser.add_argument("--token", help="API token for authentication")
    parser.add_argument("--url", help="API base URL", default="http://localhost:8000")

    args = parser.parse_args()

    if args.url:
        global API_BASE_URL, API_URL
        API_BASE_URL = args.url
        API_URL = f"{API_BASE_URL}/api/custom-scripts"

    tester = CustomScriptTester(token=args.token)
    success = tester.run_all_tests()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
