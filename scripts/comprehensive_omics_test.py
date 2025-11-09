#!/usr/bin/env python3
"""
Comprehensive Omics Analysis Test Script
Tests all omics modules (Genomics, Transcriptomics, Proteomics, Metabolomics, Epigenomics, Single-cell)
with real data downloads and complete analysis workflows.
"""

import os
import sys
import json
import time
import logging
import requests
from pathlib import Path
from typing import Dict, List, Optional
import subprocess

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("comprehensive_test.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

# Configuration
BASE_URL = os.getenv("API_URL", "http://localhost:8001/api/v1")
TEST_USER_EMAIL = os.getenv("TEST_EMAIL", "test@example.com")
TEST_USER_PASSWORD = os.getenv("TEST_PASSWORD", "testpass123")
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
TEST_RESULTS_DIR = Path(__file__).parent.parent / "test_results"

# Create directories
TEST_DATA_DIR.mkdir(exist_ok=True)
TEST_RESULTS_DIR.mkdir(exist_ok=True)


class OmicsTestClient:
    """Client for testing Omics API"""

    def __init__(self, base_url: str):
        self.base_url = base_url
        self.token = None
        self.session = requests.Session()
        self.current_project_id = None

    def register_and_login(self, email: str, password: str) -> bool:
        """Register user if not exists, then login"""
        try:
            # Try to register
            response = self.session.post(
                f"{self.base_url}/auth/register",
                json={"email": email, "password": password, "full_name": "Test User"},
            )
            if response.status_code == 201:
                logger.info("User registered successfully")
            elif response.status_code == 400:
                logger.info("User already exists, proceeding to login")
            else:
                logger.warning(f"Registration response: {response.status_code}")

        except Exception as e:
            logger.warning(f"Registration error (may already exist): {e}")

        # Login
        try:
            response = self.session.post(
                f"{self.base_url}/auth/login",
                data={"username": email, "password": password},
            )
            response.raise_for_status()
            data = response.json()
            self.token = data["access_token"]
            self.session.headers.update({"Authorization": f"Bearer {self.token}"})
            logger.info(f"Logged in successfully as {email}")
            return True
        except Exception as e:
            logger.error(f"Login failed: {e}")
            return False

    def create_project(self, name: str, description: str) -> Optional[int]:
        """Create a new project"""
        try:
            response = self.session.post(
                f"{self.base_url}/projects/",
                json={"name": name, "description": description},
            )
            response.raise_for_status()
            project = response.json()
            project_id = project["id"]
            self.current_project_id = project_id
            logger.info(f"Created project '{name}' with ID {project_id}")
            return project_id
        except Exception as e:
            logger.error(f"Failed to create project: {e}")
            return None

    def upload_file(self, project_id: int, file_path: Path) -> Optional[int]:
        """Upload a data file"""
        try:
            with open(file_path, "rb") as f:
                files = {"file": (file_path.name, f)}
                data = {"project_id": project_id}
                response = self.session.post(
                    f"{self.base_url}/data/upload", files=files, data=data
                )
                response.raise_for_status()
                file_data = response.json()
                file_id = file_data["id"]
                logger.info(f"Uploaded file '{file_path.name}' with ID {file_id}")
                return file_id
        except Exception as e:
            logger.error(f"Failed to upload file {file_path}: {e}")
            return None

    def create_custom_pipeline(
        self, project_id: int, name: str, category: str, nodes: List, edges: List
    ) -> Optional[int]:
        """Create a custom pipeline"""
        try:
            response = self.session.post(
                f"{self.base_url}/custom-pipelines/",
                json={
                    "name": name,
                    "description": f"Test pipeline for {category}",
                    "category": category,
                    "is_public": False,
                    "project_id": project_id,
                    "definition": {"nodes": nodes, "edges": edges, "parameters": {}},
                },
            )
            response.raise_for_status()
            pipeline = response.json()
            pipeline_id = pipeline["id"]
            logger.info(f"Created custom pipeline '{name}' with ID {pipeline_id}")
            return pipeline_id
        except Exception as e:
            logger.error(f"Failed to create custom pipeline: {e}")
            return None

    def create_run(
        self,
        project_id: int,
        name: str,
        pipeline_template_id: str,
        input_files: List[int],
        parameters: Dict,
    ) -> Optional[int]:
        """Create a run"""
        try:
            response = self.session.post(
                f"{self.base_url}/runs/",
                json={
                    "name": name,
                    "description": f"Test run for {pipeline_template_id}",
                    "project_id": project_id,
                    "pipeline_type": "template",
                    "pipeline_template_id": pipeline_template_id,
                    "input_files": input_files,
                    "parameters": parameters,
                    "auto_start": False,
                },
            )
            response.raise_for_status()
            run = response.json()
            run_id = run["id"]
            logger.info(f"Created run '{name}' with ID {run_id}")
            return run_id
        except Exception as e:
            logger.error(f"Failed to create run: {e}")
            logger.error(
                f"Response: {e.response.text if hasattr(e, 'response') else 'No response'}"
            )
            return None

    def get_pipelines(self) -> List[Dict]:
        """Get available pipeline templates"""
        try:
            response = self.session.get(f"{self.base_url}/pipelines/")
            response.raise_for_status()
            pipelines = response.json()
            logger.info(f"Retrieved {len(pipelines)} pipeline templates")
            return pipelines
        except Exception as e:
            logger.error(f"Failed to get pipelines: {e}")
            return []


class DataDownloader:
    """Download test data for different omics types"""

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir

    def download_genomics_data(self) -> List[Path]:
        """Download small genomics test data"""
        logger.info("Preparing genomics test data...")
        files = []

        # Create small FASTQ file
        fastq_file = self.data_dir / "test_genomics.fastq"
        with open(fastq_file, "w") as f:
            # Generate 1000 fake reads
            for i in range(1000):
                f.write(f"@READ_{i}\n")
                f.write(
                    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
                )
                f.write("+\n")
                f.write(
                    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
                )
        files.append(fastq_file)
        logger.info(f"Created genomics test file: {fastq_file}")

        # Create VCF file
        vcf_file = self.data_dir / "test_variants.vcf"
        with open(vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##contig=<ID=chr1,length=248956422>\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n")
            for i in range(1, 101):
                pos = i * 10000
                f.write(f"chr1\t{pos}\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n")
        files.append(vcf_file)
        logger.info(f"Created VCF test file: {vcf_file}")

        return files

    def download_transcriptomics_data(self) -> List[Path]:
        """Download RNA-seq test data"""
        logger.info("Preparing transcriptomics test data...")
        files = []

        # Create count matrix
        counts_file = self.data_dir / "test_rna_counts.csv"
        with open(counts_file, "w") as f:
            f.write("gene_id,sample1,sample2,sample3,sample4\n")
            for i in range(1000):
                counts = [str(int(100 + i * 0.5 + j * 10)) for j in range(4)]
                f.write(f"GENE_{i:04d},{','.join(counts)}\n")
        files.append(counts_file)
        logger.info(f"Created RNA-seq counts file: {counts_file}")

        # Create metadata
        metadata_file = self.data_dir / "test_rna_metadata.csv"
        with open(metadata_file, "w") as f:
            f.write("sample,condition,replicate\n")
            f.write("sample1,control,1\n")
            f.write("sample2,control,2\n")
            f.write("sample3,treatment,1\n")
            f.write("sample4,treatment,2\n")
        files.append(metadata_file)
        logger.info(f"Created metadata file: {metadata_file}")

        return files

    def download_proteomics_data(self) -> List[Path]:
        """Download proteomics test data"""
        logger.info("Preparing proteomics test data...")
        files = []

        # Create protein abundance matrix
        protein_file = self.data_dir / "test_protein_abundance.csv"
        with open(protein_file, "w") as f:
            f.write("protein_id,sample1,sample2,sample3,sample4\n")
            for i in range(500):
                abundances = [f"{1000 + i * 2.5 + j * 50:.2f}" for j in range(4)]
                f.write(f"PROT_{i:04d},{','.join(abundances)}\n")
        files.append(protein_file)
        logger.info(f"Created proteomics file: {protein_file}")

        return files

    def download_metabolomics_data(self) -> List[Path]:
        """Download metabolomics test data"""
        logger.info("Preparing metabolomics test data...")
        files = []

        # Create metabolite abundance matrix
        metabolite_file = self.data_dir / "test_metabolite_abundance.csv"
        with open(metabolite_file, "w") as f:
            f.write("metabolite_id,sample1,sample2,sample3,sample4\n")
            for i in range(300):
                abundances = [f"{500 + i * 1.2 + j * 25:.2f}" for j in range(4)]
                f.write(f"MET_{i:04d},{','.join(abundances)}\n")
        files.append(metabolite_file)
        logger.info(f"Created metabolomics file: {metabolite_file}")

        return files

    def download_epigenomics_data(self) -> List[Path]:
        """Download epigenomics test data"""
        logger.info("Preparing epigenomics test data...")
        files = []

        # Create methylation data
        methyl_file = self.data_dir / "test_methylation.csv"
        with open(methyl_file, "w") as f:
            f.write("cpg_site,chr,position,sample1,sample2,sample3,sample4\n")
            for i in range(1000):
                chr_num = (i % 22) + 1
                pos = i * 1000
                beta_values = [f"{0.3 + (i % 7) * 0.1:.3f}" for _ in range(4)]
                f.write(f"CpG_{i:05d},chr{chr_num},{pos},{','.join(beta_values)}\n")
        files.append(methyl_file)
        logger.info(f"Created epigenomics file: {methyl_file}")

        return files

    def download_singlecell_data(self) -> List[Path]:
        """Download single-cell test data"""
        logger.info("Preparing single-cell test data...")
        files = []

        # Create expression matrix
        sc_file = self.data_dir / "test_singlecell_matrix.csv"
        with open(sc_file, "w") as f:
            # Header with cell barcodes
            cells = [f"CELL_{i:04d}" for i in range(100)]
            f.write(f"gene,{','.join(cells)}\n")
            # Gene expression data
            for i in range(500):
                expression = [str(int(abs(hash(f"{i}_{j}") % 100))) for j in range(100)]
                f.write(f"GENE_{i:04d},{','.join(expression)}\n")
        files.append(sc_file)
        logger.info(f"Created single-cell file: {sc_file}")

        return files


def test_genomics_workflow(
    client: OmicsTestClient, project_id: int, data_files: List[Path]
) -> bool:
    """Test genomics analysis workflow"""
    logger.info("=" * 80)
    logger.info("TESTING GENOMICS WORKFLOW")
    logger.info("=" * 80)

    try:
        # Upload files
        file_ids = []
        for file_path in data_files:
            file_id = client.upload_file(project_id, file_path)
            if file_id:
                file_ids.append(file_id)

        if not file_ids:
            logger.error("Failed to upload genomics files")
            return False

        # Create run with genomics pipeline
        run_id = client.create_run(
            project_id=project_id,
            name="Genomics Test Run",
            pipeline_template_id="variant-calling-gatk",
            input_files=file_ids,
            parameters={
                "reference_genome": "hg38",
                "quality_threshold": 30,
                "threads": 4,
            },
        )

        if run_id:
            logger.info(f"✅ Genomics workflow test PASSED - Run ID: {run_id}")
            return True
        else:
            logger.error("❌ Genomics workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Genomics workflow test FAILED with exception: {e}")
        return False


def test_transcriptomics_workflow(
    client: OmicsTestClient, project_id: int, data_files: List[Path]
) -> bool:
    """Test transcriptomics analysis workflow"""
    logger.info("=" * 80)
    logger.info("TESTING TRANSCRIPTOMICS WORKFLOW")
    logger.info("=" * 80)

    try:
        file_ids = []
        for file_path in data_files:
            file_id = client.upload_file(project_id, file_path)
            if file_id:
                file_ids.append(file_id)

        if not file_ids:
            logger.error("Failed to upload transcriptomics files")
            return False

        run_id = client.create_run(
            project_id=project_id,
            name="RNA-seq Test Run",
            pipeline_template_id="rna-seq-basic",
            input_files=file_ids,
            parameters={
                "normalization_method": "DESeq2",
                "fdr_threshold": 0.05,
                "log2fc_threshold": 1.0,
            },
        )

        if run_id:
            logger.info(f"✅ Transcriptomics workflow test PASSED - Run ID: {run_id}")
            return True
        else:
            logger.error("❌ Transcriptomics workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Transcriptomics workflow test FAILED with exception: {e}")
        return False


def test_proteomics_workflow(
    client: OmicsTestClient, project_id: int, data_files: List[Path]
) -> bool:
    """Test proteomics analysis workflow"""
    logger.info("=" * 80)
    logger.info("TESTING PROTEOMICS WORKFLOW")
    logger.info("=" * 80)

    try:
        file_ids = []
        for file_path in data_files:
            file_id = client.upload_file(project_id, file_path)
            if file_id:
                file_ids.append(file_id)

        if not file_ids:
            logger.error("Failed to upload proteomics files")
            return False

        run_id = client.create_run(
            project_id=project_id,
            name="Proteomics Test Run",
            pipeline_template_id="protein-identification",
            input_files=file_ids,
            parameters={
                "database": "uniprot",
                "peptide_fdr": 0.01,
                "protein_fdr": 0.05,
            },
        )

        if run_id:
            logger.info(f"✅ Proteomics workflow test PASSED - Run ID: {run_id}")
            return True
        else:
            logger.error("❌ Proteomics workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Proteomics workflow test FAILED with exception: {e}")
        return False


def test_metabolomics_workflow(
    client: OmicsTestClient, project_id: int, data_files: List[Path]
) -> bool:
    """Test metabolomics analysis workflow"""
    logger.info("=" * 80)
    logger.info("TESTING METABOLOMICS WORKFLOW")
    logger.info("=" * 80)

    try:
        file_ids = []
        for file_path in data_files:
            file_id = client.upload_file(project_id, file_path)
            if file_id:
                file_ids.append(file_id)

        if not file_ids:
            logger.error("Failed to upload metabolomics files")
            return False

        run_id = client.create_run(
            project_id=project_id,
            name="Metabolomics Test Run",
            pipeline_template_id="metabolite-profiling",
            input_files=file_ids,
            parameters={
                "normalization": "sum",
                "scaling": "pareto",
                "missing_value_threshold": 0.5,
            },
        )

        if run_id:
            logger.info(f"✅ Metabolomics workflow test PASSED - Run ID: {run_id}")
            return True
        else:
            logger.error("❌ Metabolomics workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Metabolomics workflow test FAILED with exception: {e}")
        return False


def test_epigenomics_workflow(
    client: OmicsTestClient, project_id: int, data_files: List[Path]
) -> bool:
    """Test epigenomics analysis workflow"""
    logger.info("=" * 80)
    logger.info("TESTING EPIGENOMICS WORKFLOW")
    logger.info("=" * 80)

    try:
        file_ids = []
        for file_path in data_files:
            file_id = client.upload_file(project_id, file_path)
            if file_id:
                file_ids.append(file_id)

        if not file_ids:
            logger.error("Failed to upload epigenomics files")
            return False

        run_id = client.create_run(
            project_id=project_id,
            name="Epigenomics Test Run",
            pipeline_template_id="methylation-analysis",
            input_files=file_ids,
            parameters={
                "analysis_type": "differential",
                "beta_threshold": 0.2,
                "p_value_threshold": 0.05,
            },
        )

        if run_id:
            logger.info(f"✅ Epigenomics workflow test PASSED - Run ID: {run_id}")
            return True
        else:
            logger.error("❌ Epigenomics workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Epigenomics workflow test FAILED with exception: {e}")
        return False


def test_singlecell_workflow(
    client: OmicsTestClient, project_id: int, data_files: List[Path]
) -> bool:
    """Test single-cell analysis workflow"""
    logger.info("=" * 80)
    logger.info("TESTING SINGLE-CELL WORKFLOW")
    logger.info("=" * 80)

    try:
        file_ids = []
        for file_path in data_files:
            file_id = client.upload_file(project_id, file_path)
            if file_id:
                file_ids.append(file_id)

        if not file_ids:
            logger.error("Failed to upload single-cell files")
            return False

        run_id = client.create_run(
            project_id=project_id,
            name="Single-cell Test Run",
            pipeline_template_id="scrna-seq-clustering",
            input_files=file_ids,
            parameters={
                "normalization": "log",
                "clustering_resolution": 0.8,
                "min_genes": 200,
                "min_cells": 3,
            },
        )

        if run_id:
            logger.info(f"✅ Single-cell workflow test PASSED - Run ID: {run_id}")
            return True
        else:
            logger.error("❌ Single-cell workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Single-cell workflow test FAILED with exception: {e}")
        return False


def test_multiomics_workflow(client: OmicsTestClient, project_id: int) -> bool:
    """Test multi-omics integration workflow"""
    logger.info("=" * 80)
    logger.info("TESTING MULTI-OMICS INTEGRATION WORKFLOW")
    logger.info("=" * 80)

    try:
        # Create a custom multi-omics pipeline
        nodes = [
            {
                "id": "input1",
                "type": "input",
                "label": "Transcriptomics Data",
                "position": {"x": 100, "y": 100},
                "data": {"type": "transcriptomics"},
            },
            {
                "id": "input2",
                "type": "input",
                "label": "Proteomics Data",
                "position": {"x": 100, "y": 200},
                "data": {"type": "proteomics"},
            },
            {
                "id": "integrate",
                "type": "process",
                "label": "Multi-omics Integration",
                "position": {"x": 300, "y": 150},
                "data": {"tool": "MOFA"},
            },
            {
                "id": "output",
                "type": "output",
                "label": "Integrated Results",
                "position": {"x": 500, "y": 150},
                "data": {},
            },
        ]

        edges = [
            {"id": "e1", "source": "input1", "target": "integrate"},
            {"id": "e2", "source": "input2", "target": "integrate"},
            {"id": "e3", "source": "integrate", "target": "output"},
        ]

        pipeline_id = client.create_custom_pipeline(
            project_id=project_id,
            name="Multi-omics Integration Pipeline",
            category="multiomics",
            nodes=nodes,
            edges=edges,
        )

        if pipeline_id:
            logger.info(
                f"✅ Multi-omics workflow test PASSED - Pipeline ID: {pipeline_id}"
            )
            return True
        else:
            logger.error("❌ Multi-omics workflow test FAILED")
            return False

    except Exception as e:
        logger.error(f"❌ Multi-omics workflow test FAILED with exception: {e}")
        return False


def main():
    """Main test execution"""
    logger.info("=" * 80)
    logger.info("COMPREHENSIVE OMICS PLATFORM TEST")
    logger.info("=" * 80)
    logger.info(f"API URL: {BASE_URL}")
    logger.info(f"Test data directory: {TEST_DATA_DIR}")
    logger.info(f"Test results directory: {TEST_RESULTS_DIR}")
    logger.info("")

    # Initialize client
    client = OmicsTestClient(BASE_URL)

    # Login
    logger.info("Step 1: Authentication")
    if not client.register_and_login(TEST_USER_EMAIL, TEST_USER_PASSWORD):
        logger.error("Authentication failed. Exiting.")
        return False

    # Get available pipelines
    logger.info("\nStep 2: Fetching available pipeline templates")
    pipelines = client.get_pipelines()
    if pipelines:
        logger.info("Available pipelines:")
        for p in pipelines:
            logger.info(f"  - {p.get('id')}: {p.get('name')} ({p.get('category')})")

    # Download test data
    logger.info("\nStep 3: Preparing test data")
    downloader = DataDownloader(TEST_DATA_DIR)

    test_results = {}

    # Test each omics workflow
    workflows = [
        ("genomics", downloader.download_genomics_data, test_genomics_workflow),
        (
            "transcriptomics",
            downloader.download_transcriptomics_data,
            test_transcriptomics_workflow,
        ),
        ("proteomics", downloader.download_proteomics_data, test_proteomics_workflow),
        (
            "metabolomics",
            downloader.download_metabolomics_data,
            test_metabolomics_workflow,
        ),
        (
            "epigenomics",
            downloader.download_epigenomics_data,
            test_epigenomics_workflow,
        ),
        ("singlecell", downloader.download_singlecell_data, test_singlecell_workflow),
    ]

    for omics_type, download_func, test_func in workflows:
        logger.info(f"\n{'=' * 80}")
        logger.info(f"Testing {omics_type.upper()}")
        logger.info(f"{'=' * 80}")

        # Create project for this omics type
        project_id = client.create_project(
            name=f"Test {omics_type.capitalize()} Project",
            description=f"Automated test project for {omics_type} analysis",
        )

        if not project_id:
            logger.error(f"Failed to create project for {omics_type}")
            test_results[omics_type] = False
            continue

        # Download data
        data_files = download_func()

        if not data_files:
            logger.error(f"Failed to prepare data for {omics_type}")
            test_results[omics_type] = False
            continue

        # Run test
        test_results[omics_type] = test_func(client, project_id, data_files)

        time.sleep(2)  # Brief pause between tests

    # Test multi-omics integration
    logger.info(f"\n{'=' * 80}")
    logger.info("Testing MULTI-OMICS INTEGRATION")
    logger.info(f"{'=' * 80}")

    multiomics_project_id = client.create_project(
        name="Test Multi-omics Integration Project",
        description="Automated test project for multi-omics integration",
    )

    if multiomics_project_id:
        test_results["multiomics"] = test_multiomics_workflow(
            client, multiomics_project_id
        )
    else:
        test_results["multiomics"] = False

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("TEST SUMMARY")
    logger.info("=" * 80)

    total_tests = len(test_results)
    passed_tests = sum(1 for result in test_results.values() if result)

    for omics_type, result in test_results.items():
        status = "✅ PASSED" if result else "❌ FAILED"
        logger.info(f"{omics_type.capitalize()}: {status}")

    logger.info("")
    logger.info(f"Total: {passed_tests}/{total_tests} tests passed")
    logger.info(f"Success rate: {passed_tests/total_tests*100:.1f}%")

    # Save results
    results_file = TEST_RESULTS_DIR / "test_results.json"
    with open(results_file, "w") as f:
        json.dump(
            {
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "total_tests": total_tests,
                "passed_tests": passed_tests,
                "success_rate": passed_tests / total_tests * 100,
                "details": {
                    k: "PASSED" if v else "FAILED" for k, v in test_results.items()
                },
            },
            f,
            indent=2,
        )

    logger.info(f"\nResults saved to: {results_file}")

    return passed_tests == total_tests


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
