#!/usr/bin/env python3
"""
Comprehensive test script for Visualization Workspace functionality.

Tests:
1. Dashboard CRUD operations
2. Panel management
3. Layout system
4. Template system
5. Batch operations
6. Error handling
"""

import sys
import json
import requests
from typing import Optional, Dict, Any, List

# Configuration
BASE_URL = "http://localhost:8000"
API_BASE = f"{BASE_URL}/api/viz-workspace"

# Test data
TEST_PROJECT_ID = 1


class Colors:
    """ANSI color codes for terminal output"""

    GREEN = "\033[92m"
    RED = "\033[91m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    RESET = "\033[0m"
    BOLD = "\033[1m"


def print_section(title: str):
    """Print a section header"""
    print(f"\n{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.RESET}")
    print(f"{Colors.CYAN}{Colors.BOLD}{title}{Colors.RESET}")
    print(f"{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.RESET}\n")


def print_test(test_name: str):
    """Print test name"""
    print(f"{Colors.BLUE}Testing: {test_name}{Colors.RESET}")


def print_success(message: str):
    """Print success message"""
    print(f"{Colors.GREEN}✓ {message}{Colors.RESET}")


def print_error(message: str):
    """Print error message"""
    print(f"{Colors.RED}✗ {message}{Colors.RESET}")


def print_warning(message: str):
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠ {message}{Colors.RESET}")


def print_info(message: str):
    """Print info message"""
    print(f"{Colors.CYAN}ℹ {message}{Colors.RESET}")


def print_json(data: Dict[Any, Any], indent: int = 2):
    """Pretty print JSON data"""
    print(json.dumps(data, indent=indent, default=str))


class VisualizationWorkspaceTest:
    """Test suite for Visualization Workspace"""

    def __init__(self):
        self.session = requests.Session()
        self.test_dashboard_id: Optional[int] = None
        self.test_panel_ids: List[int] = []
        self.test_template_id: Optional[int] = None
        self.passed = 0
        self.failed = 0

    def assert_response(
        self, response: requests.Response, expected_status: int, test_name: str
    ) -> bool:
        """Assert response status and handle result"""
        if response.status_code == expected_status:
            self.passed += 1
            print_success(f"{test_name} - Status {response.status_code}")
            return True
        else:
            self.failed += 1
            print_error(
                f"{test_name} - Expected {expected_status}, got {response.status_code}"
            )
            try:
                error_data = response.json()
                print_info(f"Error: {error_data}")
            except:
                print_info(f"Response: {response.text}")
            return False

    def test_create_dashboard(self) -> bool:
        """Test: Create a new dashboard"""
        print_test("Create Dashboard")

        data = {
            "project_id": TEST_PROJECT_ID,
            "name": "Test RNA-seq Dashboard",
            "description": "Test dashboard for RNA-seq analysis",
            "layout": {
                "cols": 12,
                "rowHeight": 100,
                "breakpoints": {"lg": 1200, "md": 996, "sm": 768},
            },
            "metadata": {"experiment": "RNA-seq", "version": "1.0"},
            "is_public": False,
        }

        response = self.session.post(f"{API_BASE}/dashboards", json=data)

        if self.assert_response(response, 200, "Create dashboard"):
            result = response.json()
            self.test_dashboard_id = result["id"]
            print_info(f"Created dashboard ID: {self.test_dashboard_id}")
            print_info(f"Dashboard name: {result['name']}")
            return True
        return False

    def test_get_dashboard(self) -> bool:
        """Test: Get dashboard by ID"""
        print_test("Get Dashboard")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        response = self.session.get(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}",
            params={"include_panels": True},
        )

        if self.assert_response(response, 200, "Get dashboard"):
            result = response.json()
            print_info(f"Dashboard: {result['name']}")
            print_info(f"Panels: {result.get('panel_count', 0)}")
            return True
        return False

    def test_list_dashboards(self) -> bool:
        """Test: List all dashboards"""
        print_test("List Dashboards")

        response = self.session.get(
            f"{API_BASE}/dashboards", params={"project_id": TEST_PROJECT_ID}
        )

        if self.assert_response(response, 200, "List dashboards"):
            result = response.json()
            print_info(f"Found {len(result)} dashboard(s)")
            for dashboard in result[:3]:  # Show first 3
                print_info(f"  - {dashboard['name']} (ID: {dashboard['id']})")
            return True
        return False

    def test_add_panel_heatmap(self) -> bool:
        """Test: Add heatmap panel"""
        print_test("Add Heatmap Panel")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {
            "panel_key": "expression_heatmap",
            "title": "Gene Expression Heatmap",
            "description": "Top 50 differentially expressed genes",
            "viz_type": "heatmap",
            "data_source": {
                "type": "dataset",
                "dataset_id": 123,
                "file_path": "/data/expression_matrix.csv",
            },
            "viz_config": {
                "colormap": "RdBu",
                "normalize": "row",
                "cluster_rows": True,
                "show_legend": True,
            },
            "position": {"x": 0, "y": 0, "w": 6, "h": 4},
            "auto_refresh": False,
            "refresh_interval": None,
        }

        response = self.session.post(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels", json=data
        )

        if self.assert_response(response, 200, "Add heatmap panel"):
            result = response.json()
            self.test_panel_ids.append(result["id"])
            print_info(f"Created panel ID: {result['id']}")
            print_info(f"Panel key: {result['panel_key']}")
            return True
        return False

    def test_add_panel_line_chart(self) -> bool:
        """Test: Add line chart panel"""
        print_test("Add Line Chart Panel")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {
            "panel_key": "expression_trend",
            "title": "Expression Trend",
            "description": "Gene expression over time",
            "viz_type": "line",
            "data_source": {"type": "file", "file_path": "/data/time_series.csv"},
            "viz_config": {
                "x_column": "time",
                "y_columns": ["gene1", "gene2"],
                "line_style": "solid",
                "grid": True,
            },
            "position": {"x": 6, "y": 0, "w": 6, "h": 4},
            "auto_refresh": True,
            "refresh_interval": 60,
        }

        response = self.session.post(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels", json=data
        )

        if self.assert_response(response, 200, "Add line chart panel"):
            result = response.json()
            self.test_panel_ids.append(result["id"])
            print_info(f"Created panel ID: {result['id']}")
            print_info(f"Auto-refresh: {result['auto_refresh']}")
            return True
        return False

    def test_add_panel_metric(self) -> bool:
        """Test: Add metric panel"""
        print_test("Add Metric Panel")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {
            "panel_key": "sample_count",
            "title": "Sample Count",
            "description": "Total number of samples",
            "viz_type": "metric",
            "data_source": {"type": "api", "url": "/api/samples/count"},
            "viz_config": {"value_format": ".0f", "color": "#4CAF50"},
            "position": {"x": 0, "y": 4, "w": 3, "h": 2},
        }

        response = self.session.post(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels", json=data
        )

        if self.assert_response(response, 200, "Add metric panel"):
            result = response.json()
            self.test_panel_ids.append(result["id"])
            print_info(f"Created panel ID: {result['id']}")
            return True
        return False

    def test_list_panels(self) -> bool:
        """Test: List all panels in dashboard"""
        print_test("List Dashboard Panels")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        response = self.session.get(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels"
        )

        if self.assert_response(response, 200, "List panels"):
            result = response.json()
            print_info(f"Found {len(result)} panel(s)")
            for panel in result:
                print_info(f"  - {panel['title']} ({panel['viz_type']})")
            return True
        return False

    def test_update_panel(self) -> bool:
        """Test: Update panel configuration"""
        print_test("Update Panel")

        if not self.test_panel_ids:
            print_warning("No test panels available, skipping")
            return False

        panel_id = self.test_panel_ids[0]
        data = {
            "title": "Updated Gene Expression Heatmap",
            "description": "Updated description",
            "viz_config": {
                "colormap": "viridis",
                "normalize": "column",
                "cluster_rows": False,
                "show_legend": True,
            },
            "position": {"x": 0, "y": 0, "w": 8, "h": 5},
        }

        response = self.session.put(f"{API_BASE}/panels/{panel_id}", json=data)

        if self.assert_response(response, 200, "Update panel"):
            result = response.json()
            print_info(f"Updated panel: {result['title']}")
            print_info(
                f"New position: w={result['position']['w']}, h={result['position']['h']}"
            )
            return True
        return False

    def test_batch_update_positions(self) -> bool:
        """Test: Batch update panel positions"""
        print_test("Batch Update Panel Positions")

        if not self.test_dashboard_id or len(self.test_panel_ids) < 2:
            print_warning("Insufficient test data, skipping")
            return False

        # Get panels to get their keys
        panels_response = self.session.get(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels"
        )

        if panels_response.status_code != 200:
            print_error("Failed to get panels for batch update")
            return False

        panels = panels_response.json()
        positions = {}

        # Rearrange panels
        for i, panel in enumerate(panels):
            positions[panel["panel_key"]] = {
                "x": (i * 4) % 12,
                "y": (i // 3) * 4,
                "w": 4,
                "h": 3,
            }

        data = {"positions": positions}

        response = self.session.put(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels/positions",
            json=data,
        )

        if self.assert_response(response, 200, "Batch update positions"):
            print_info(f"Updated {len(positions)} panel positions")
            return True
        return False

    def test_dashboard_statistics(self) -> bool:
        """Test: Get dashboard statistics"""
        print_test("Get Dashboard Statistics")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        response = self.session.get(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/stats"
        )

        if self.assert_response(response, 200, "Get statistics"):
            result = response.json()
            print_info(f"Panel count: {result['panel_count']}")
            print_info(f"Visualization types: {result['viz_types']}")
            print_info(f"Auto-refresh panels: {result['auto_refresh_panels']}")
            return True
        return False

    def test_duplicate_dashboard(self) -> bool:
        """Test: Duplicate dashboard"""
        print_test("Duplicate Dashboard")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {"new_name": "Copy of Test RNA-seq Dashboard"}

        response = self.session.post(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/duplicate", json=data
        )

        if self.assert_response(response, 200, "Duplicate dashboard"):
            result = response.json()
            print_info(f"Duplicated dashboard ID: {result['id']}")
            print_info(f"New name: {result['name']}")
            print_info(f"Copied {result.get('panel_count', 0)} panels")

            # Clean up duplicate
            cleanup_response = self.session.delete(
                f"{API_BASE}/dashboards/{result['id']}"
            )
            if cleanup_response.status_code == 204:
                print_info("Cleaned up duplicate dashboard")

            return True
        return False

    def test_update_dashboard(self) -> bool:
        """Test: Update dashboard properties"""
        print_test("Update Dashboard")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {
            "name": "Updated RNA-seq Dashboard",
            "description": "Updated description for testing",
            "is_public": True,
            "metadata": {"experiment": "RNA-seq", "version": "2.0", "updated": True},
        }

        response = self.session.put(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}", json=data
        )

        if self.assert_response(response, 200, "Update dashboard"):
            result = response.json()
            print_info(f"Updated dashboard: {result['name']}")
            print_info(f"Is public: {result['is_public']}")
            return True
        return False

    def test_create_template(self) -> bool:
        """Test: Create template from dashboard"""
        print_test("Create Template from Dashboard")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {"template_name": "Test RNA-seq Template"}

        response = self.session.post(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/make-template", json=data
        )

        if self.assert_response(response, 200, "Create template"):
            result = response.json()
            self.test_template_id = result["id"]
            print_info(f"Created template ID: {self.test_template_id}")
            print_info(f"Template name: {result['name']}")
            print_info(f"Is template: {result.get('is_template', False)}")
            return True
        return False

    def test_list_templates(self) -> bool:
        """Test: List all templates"""
        print_test("List Templates")

        response = self.session.get(f"{API_BASE}/templates")

        if self.assert_response(response, 200, "List templates"):
            result = response.json()
            print_info(f"Found {len(result)} template(s)")
            for template in result[:3]:  # Show first 3
                print_info(f"  - {template['name']} (ID: {template['id']})")
            return True
        return False

    def test_create_from_template(self) -> bool:
        """Test: Create dashboard from template"""
        print_test("Create Dashboard from Template")

        if not self.test_template_id:
            print_warning("No test template available, skipping")
            return False

        data = {
            "template_id": self.test_template_id,
            "name": "New Dashboard from Template",
            "project_id": TEST_PROJECT_ID,
        }

        response = self.session.post(
            f"{API_BASE}/templates/create-dashboard", json=data
        )

        if self.assert_response(response, 200, "Create from template"):
            result = response.json()
            print_info(f"Created dashboard ID: {result['id']}")
            print_info(f"Dashboard name: {result['name']}")
            print_info(f"Panels copied: {result.get('panel_count', 0)}")

            # Clean up created dashboard
            cleanup_response = self.session.delete(
                f"{API_BASE}/dashboards/{result['id']}"
            )
            if cleanup_response.status_code == 204:
                print_info("Cleaned up template-created dashboard")

            return True
        return False

    def test_delete_panel(self) -> bool:
        """Test: Delete a panel"""
        print_test("Delete Panel")

        if not self.test_panel_ids:
            print_warning("No test panels available, skipping")
            return False

        panel_id = self.test_panel_ids[-1]  # Delete last panel

        response = self.session.delete(f"{API_BASE}/panels/{panel_id}")

        if self.assert_response(response, 204, "Delete panel"):
            print_info(f"Deleted panel ID: {panel_id}")
            self.test_panel_ids.pop()
            return True
        return False

    def test_error_handling_invalid_dashboard(self) -> bool:
        """Test: Error handling for invalid dashboard ID"""
        print_test("Error Handling - Invalid Dashboard ID")

        response = self.session.get(f"{API_BASE}/dashboards/999999")

        if self.assert_response(response, 404, "Invalid dashboard ID"):
            print_info("Correctly returned 404 for non-existent dashboard")
            return True
        return False

    def test_error_handling_duplicate_panel_key(self) -> bool:
        """Test: Error handling for duplicate panel key"""
        print_test("Error Handling - Duplicate Panel Key")

        if not self.test_dashboard_id:
            print_warning("No test dashboard available, skipping")
            return False

        data = {
            "panel_key": "expression_heatmap",  # Already exists
            "title": "Duplicate Panel",
            "viz_type": "line",
            "position": {"x": 0, "y": 0, "w": 4, "h": 3},
        }

        response = self.session.post(
            f"{API_BASE}/dashboards/{self.test_dashboard_id}/panels", json=data
        )

        if self.assert_response(response, 400, "Duplicate panel key"):
            print_info("Correctly rejected duplicate panel key")
            return True
        return False

    def cleanup(self) -> bool:
        """Clean up test data"""
        print_section("Cleanup Test Data")

        success = True

        # Delete template
        if self.test_template_id:
            response = self.session.delete(
                f"{API_BASE}/dashboards/{self.test_template_id}"
            )
            if response.status_code == 204:
                print_success(f"Deleted template ID: {self.test_template_id}")
            else:
                print_warning(f"Failed to delete template: {response.status_code}")
                success = False

        # Delete dashboard (cascades to panels)
        if self.test_dashboard_id:
            response = self.session.delete(
                f"{API_BASE}/dashboards/{self.test_dashboard_id}"
            )
            if response.status_code == 204:
                print_success(f"Deleted dashboard ID: {self.test_dashboard_id}")
            else:
                print_warning(f"Failed to delete dashboard: {response.status_code}")
                success = False

        return success

    def run_all_tests(self):
        """Run all tests in sequence"""
        print_section("Visualization Workspace Test Suite")

        tests = [
            # Dashboard tests
            self.test_create_dashboard,
            self.test_get_dashboard,
            self.test_list_dashboards,
            # Panel tests
            self.test_add_panel_heatmap,
            self.test_add_panel_line_chart,
            self.test_add_panel_metric,
            self.test_list_panels,
            self.test_update_panel,
            self.test_batch_update_positions,
            # Statistics
            self.test_dashboard_statistics,
            # Dashboard operations
            self.test_update_dashboard,
            self.test_duplicate_dashboard,
            # Template tests
            self.test_create_template,
            self.test_list_templates,
            self.test_create_from_template,
            # Delete tests
            self.test_delete_panel,
            # Error handling
            self.test_error_handling_invalid_dashboard,
            self.test_error_handling_duplicate_panel_key,
        ]

        for test in tests:
            try:
                test()
            except Exception as e:
                self.failed += 1
                print_error(f"Test failed with exception: {str(e)}")
            print()  # Blank line between tests

        # Cleanup
        self.cleanup()

        # Print summary
        self.print_summary()

    def print_summary(self):
        """Print test summary"""
        print_section("Test Summary")

        total = self.passed + self.failed
        success_rate = (self.passed / total * 100) if total > 0 else 0

        print(f"Total tests: {total}")
        print(f"{Colors.GREEN}Passed: {self.passed}{Colors.RESET}")
        print(f"{Colors.RED}Failed: {self.failed}{Colors.RESET}")
        print(f"Success rate: {success_rate:.1f}%")

        if self.failed == 0:
            print(f"\n{Colors.GREEN}{Colors.BOLD}✓ All tests passed!{Colors.RESET}\n")
        else:
            print(f"\n{Colors.RED}{Colors.BOLD}✗ Some tests failed{Colors.RESET}\n")


def main():
    """Main test runner"""
    print(f"{Colors.BOLD}Visualization Workspace Test Suite{Colors.RESET}")
    print(f"Testing against: {BASE_URL}\n")

    # Check if API is available
    try:
        response = requests.get(f"{BASE_URL}/docs")
        if response.status_code != 200:
            print_error(f"API not available at {BASE_URL}")
            print_info("Please ensure the backend is running")
            sys.exit(1)
    except requests.exceptions.ConnectionError:
        print_error(f"Cannot connect to {BASE_URL}")
        print_info("Please ensure the backend is running")
        sys.exit(1)

    # Run tests
    test_suite = VisualizationWorkspaceTest()
    test_suite.run_all_tests()

    # Exit with appropriate code
    sys.exit(0 if test_suite.failed == 0 else 1)


if __name__ == "__main__":
    main()
