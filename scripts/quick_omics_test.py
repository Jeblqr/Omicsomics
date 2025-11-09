#!/usr/bin/env python3
"""
Quick Omics Platform Test
Tests basic functionality of all modules with correct API endpoints.
"""

import requests
import json
import time
from pathlib import Path
import sys

# Configuration
BASE_URL = "http://localhost:8001/api/v1"
TEST_USER = {
    "email": "test_user@omics.com",
    "password": "TestPassword123!",
    "full_name": "Test User",
}

# Directories
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
DOWNLOADS_DIR = PROJECT_ROOT / "downloads"
TEST_RESULTS_DIR = PROJECT_ROOT / "test_results"

DOWNLOADS_DIR.mkdir(exist_ok=True)
TEST_RESULTS_DIR.mkdir(exist_ok=True)


class Colors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"


def print_header(text):
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}\n")


def print_success(text):
    print(f"{Colors.OKGREEN}✅ {text}{Colors.ENDC}")


def print_error(text):
    print(f"{Colors.FAIL}❌ {text}{Colors.ENDC}")


def print_info(text):
    print(f"{Colors.OKBLUE}ℹ️  {text}{Colors.ENDC}")


def print_warning(text):
    print(f"{Colors.WARNING}⚠️  {text}{Colors.ENDC}")


class OmicsAPI:
    def __init__(self):
        self.base_url = BASE_URL
        self.token = None
        self.headers = {}

    def register(self):
        """Register new user"""
        print_info(f"Registering user: {TEST_USER['email']}")
        try:
            response = requests.post(
                f"{self.base_url}/auth/register",
                json={
                    "email": TEST_USER["email"],
                    "password": TEST_USER["password"],
                    "full_name": TEST_USER["full_name"],
                },
            )
            if response.status_code == 201:
                print_success("User registered successfully")
                return response.json()
            elif response.status_code == 400:
                print_warning("User already exists")
                return None
            else:
                print_error(f"Registration failed: {response.status_code}")
                print_error(response.text)
                return None
        except Exception as e:
            print_error(f"Registration error: {e}")
            return None

    def login(self):
        """Login and get token"""
        print_info(f"Logging in as: {TEST_USER['email']}")
        try:
            response = requests.post(
                f"{self.base_url}/auth/login/access-token",
                data={
                    "username": TEST_USER["email"],
                    "password": TEST_USER["password"],
                },
            )
            response.raise_for_status()
            data = response.json()
            self.token = data["access_token"]
            self.headers = {"Authorization": f"Bearer {self.token}"}
            print_success("Login successful")
            return True
        except Exception as e:
            print_error(f"Login failed: {e}")
            return False

    def create_project(self, name, description):
        """Create a project"""
        print_info(f"Creating project: {name}")
        try:
            response = requests.post(
                f"{self.base_url}/projects/",
                json={"name": name, "description": description},
                headers=self.headers,
            )
            response.raise_for_status()
            project = response.json()
            print_success(f"Project created with ID: {project['id']}")
            return project
        except Exception as e:
            print_error(f"Failed to create project: {e}")
            return None

    def upload_file(self, project_id, file_path):
        """Upload a file"""
        print_info(f"Uploading file: {file_path.name}")
        try:
            with open(file_path, "rb") as f:
                files = {"file": (file_path.name, f)}
                data = {"project_id": str(project_id)}
                response = requests.post(
                    f"{self.base_url}/data/upload",
                    files=files,
                    data=data,
                    headers=self.headers,
                )
            response.raise_for_status()
            file_data = response.json()
            print_success(f"File uploaded with ID: {file_data['id']}")
            return file_data
        except Exception as e:
            print_error(f"Failed to upload file: {e}")
            return None

    def list_pipelines(self):
        """List pipeline templates"""
        try:
            response = requests.get(f"{self.base_url}/pipelines/", headers=self.headers)
            response.raise_for_status()
            pipelines = response.json()
            print_success(f"Found {len(pipelines)} pipeline templates")
            return pipelines
        except Exception as e:
            print_error(f"Failed to list pipelines: {e}")
            return []

    def create_run(
        self,
        project_id,
        name,
        pipeline_type,
        template_id=None,
        custom_id=None,
        input_files=None,
        parameters=None,
    ):
        """Create a run"""
        print_info(f"Creating run: {name}")
        try:
            payload = {
                "name": name,
                "description": f"Test run for {name}",
                "project_id": project_id,
                "pipeline_type": pipeline_type,
                "pipeline_template_id": template_id,
                "custom_pipeline_id": custom_id,
                "parameters": parameters or {},
                "input_files": input_files or [],
                "input_mapping": {},
                "auto_start": False,
                "priority": 0,
            }
            response = requests.post(
                f"{self.base_url}/runs/", json=payload, headers=self.headers
            )
            response.raise_for_status()
            run = response.json()
            print_success(f"Run created with ID: {run['id']}")
            return run
        except Exception as e:
            print_error(f"Failed to create run: {e}")
            if hasattr(e, "response"):
                print_error(f"Response: {e.response.text}")
            return None

    def start_run(self, run_id):
        """Start a run"""
        print_info(f"Starting run ID: {run_id}")
        try:
            response = requests.post(
                f"{self.base_url}/runs/{run_id}/start", headers=self.headers
            )
            response.raise_for_status()
            print_success("Run started")
            return response.json()
        except Exception as e:
            print_error(f"Failed to start run: {e}")
            return None


def create_test_file(filename, content):
    """Create a test file"""
    file_path = DOWNLOADS_DIR / filename
    with open(file_path, "w") as f:
        f.write(content)
    print_success(f"Created test file: {filename}")
    return file_path


def test_genomics(api):
    """Test genomics workflow"""
    print_header("🧬 TESTING GENOMICS")

    # Create project
    project = api.create_project("Genomics Test", "Testing variant calling pipeline")
    if not project:
        return False

    # Create test VCF file
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	.
chr1	200	.	C	T	40	PASS	.
chr1	300	.	G	A	50	PASS	.
"""
    vcf_file = create_test_file("test_variants.vcf", vcf_content)

    # Upload file
    uploaded_file = api.upload_file(project["id"], vcf_file)
    if not uploaded_file:
        return False

    # List pipelines
    pipelines = api.list_pipelines()
    variant_pipeline = next(
        (p for p in pipelines if "variant" in p["id"].lower()), None
    )

    if variant_pipeline:
        # Create run
        run = api.create_run(
            project["id"],
            "Variant Calling Test",
            "template",
            template_id=variant_pipeline["id"],
            input_files=[uploaded_file["id"]],
            parameters={"min_quality": "30"},
        )

        if run:
            # Start run
            api.start_run(run["id"])
            print_success("Genomics test completed")
            return True
    else:
        print_warning("No variant calling pipeline found")

    return False


def test_transcriptomics(api):
    """Test transcriptomics workflow"""
    print_header("📊 TESTING TRANSCRIPTOMICS")

    # Create project
    project = api.create_project(
        "RNA-seq Test", "Testing differential expression pipeline"
    )
    if not project:
        return False

    # Create test counts file
    counts_content = """Gene	Sample1	Sample2	Sample3	Sample4
Gene1	100	110	95	105
Gene2	200	210	195	205
Gene3	50	55	48	52
Gene4	150	160	145	155
Gene5	300	320	290	310
"""
    counts_file = create_test_file("test_counts.tsv", counts_content)

    # Upload file
    uploaded_file = api.upload_file(project["id"], counts_file)
    if not uploaded_file:
        return False

    # List pipelines
    pipelines = api.list_pipelines()
    rna_pipeline = next((p for p in pipelines if "rna" in p["id"].lower()), None)

    if rna_pipeline:
        # Create run
        run = api.create_run(
            project["id"],
            "DEG Analysis Test",
            "template",
            template_id=rna_pipeline["id"],
            input_files=[uploaded_file["id"]],
            parameters={"fdr_threshold": "0.05"},
        )

        if run:
            api.start_run(run["id"])
            print_success("Transcriptomics test completed")
            return True
    else:
        print_warning("No RNA-seq pipeline found")

    return False


def test_proteomics(api):
    """Test proteomics workflow"""
    print_header("🧪 TESTING PROTEOMICS")

    # Create project
    project = api.create_project(
        "Proteomics Test", "Testing protein quantification pipeline"
    )
    if not project:
        return False

    # Create test protein data
    protein_content = """Protein,Sample1,Sample2,Sample3,Sample4
PROT_001,100,110,95,105
PROT_002,200,210,195,205
PROT_003,50,55,48,52
PROT_004,150,160,145,155
PROT_005,300,320,290,310
"""
    protein_file = create_test_file("test_proteins.csv", protein_content)

    # Upload file
    uploaded_file = api.upload_file(project["id"], protein_file)
    if not uploaded_file:
        return False

    # List pipelines
    pipelines = api.list_pipelines()
    protein_pipeline = next(
        (p for p in pipelines if "protein" in p["id"].lower()), None
    )

    if protein_pipeline:
        run = api.create_run(
            project["id"],
            "Protein Quantification Test",
            "template",
            template_id=protein_pipeline["id"],
            input_files=[uploaded_file["id"]],
            parameters={"normalization": "median"},
        )

        if run:
            api.start_run(run["id"])
            print_success("Proteomics test completed")
            return True
    else:
        print_warning("No proteomics pipeline found")

    return False


def main():
    """Main test execution"""
    print_header("🧬 OMICSOMICS PLATFORM TEST SUITE 🧬")

    api = OmicsAPI()

    # Authentication
    print_header("🔐 AUTHENTICATION")
    api.register()
    if not api.login():
        print_error("Authentication failed. Cannot continue.")
        return 1

    # Test each omics type
    results = {
        "genomics": test_genomics(api),
        "transcriptomics": test_transcriptomics(api),
        "proteomics": test_proteomics(api),
    }

    # Summary
    print_header("📋 TEST SUMMARY")
    total = len(results)
    passed = sum(results.values())

    for test_name, passed_test in results.items():
        status = "PASSED" if passed_test else "FAILED"
        color = Colors.OKGREEN if passed_test else Colors.FAIL
        print(f"{color}{test_name.upper()}: {status}{Colors.ENDC}")

    print(f"\n{Colors.BOLD}Total: {passed}/{total} tests passed{Colors.ENDC}\n")

    # Save results
    result_file = TEST_RESULTS_DIR / f"test_results_{int(time.time())}.json"
    with open(result_file, "w") as f:
        json.dump(
            {
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "results": results,
                "passed": passed,
                "total": total,
            },
            f,
            indent=2,
        )

    print_info(f"Results saved to: {result_file}")

    if passed == total:
        print_success("ALL TESTS PASSED!")
        return 0
    else:
        print_error(f"{total - passed} TEST(S) FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
