"""
Single-Cell Data Import Wizard

This scenario provides an interactive wizard for importing single-cell RNA-seq data
from various formats (10X Genomics, h5ad, Loom) and creating analysis-ready objects
(Seurat RDS, AnnData h5ad).

Key Features:
- File format detection (10X MTX, 10X HDF5, h5ad, Loom, CSV)
- Species and genome version selection
- QC threshold configuration (min_cells, min_genes, max_mito%)
- Filtering preview with statistics
- Output format selection (Seurat RDS, AnnData h5ad)
- Metadata integration
"""

from typing import Dict, List, Optional, Any, Callable
import pandas as pd
import numpy as np
from pathlib import Path
from enum import Enum
import json

from ..interactive_converter import (
    ConversionScenario,
    ConversionParameter,
    ParameterType,
    ValidationMessage,
    ValidationLevel,
    ConversionPreview,
    ConversionProgress,
)


class SingleCellInputFormat(str, Enum):
    """Single-cell input formats"""

    MTX_10X = "mtx_10x"  # 10X Genomics MTX format (matrix.mtx, genes.tsv, barcodes.tsv)
    HDF5_10X = "hdf5_10x"  # 10X Genomics HDF5 format (filtered_feature_bc_matrix.h5)
    H5AD = "h5ad"  # AnnData HDF5 format
    LOOM = "loom"  # Loom format
    CSV = "csv"  # CSV matrix (genes x cells or cells x genes)


class Species(str, Enum):
    """Supported species"""

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"
    ZEBRAFISH = "zebrafish"
    FLY = "fly"
    WORM = "worm"


class GenomeVersion(str, Enum):
    """Genome versions"""

    # Human
    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    HG38 = "hg38"
    HG19 = "hg19"

    # Mouse
    GRCM39 = "GRCm39"
    GRCM38 = "GRCm38"
    MM10 = "mm10"
    MM9 = "mm9"

    # Rat
    RNOR6 = "Rnor_6.0"
    RN6 = "rn6"


class OutputFormat(str, Enum):
    """Output format options"""

    H5AD = "h5ad"  # AnnData (Python/scanpy)
    RDS = "rds"  # Seurat (R)
    LOOM = "loom"  # Loom format


class SingleCellImportScenario(ConversionScenario):
    """Single-cell data import wizard scenario"""

    # Genome version mappings
    GENOME_VERSIONS = {
        Species.HUMAN: [
            GenomeVersion.GRCH38,
            GenomeVersion.GRCH37,
            GenomeVersion.HG38,
            GenomeVersion.HG19,
        ],
        Species.MOUSE: [
            GenomeVersion.GRCM39,
            GenomeVersion.GRCM38,
            GenomeVersion.MM10,
            GenomeVersion.MM9,
        ],
        Species.RAT: [GenomeVersion.RNOR6, GenomeVersion.RN6],
    }

    def __init__(self):
        super().__init__(
            scenario_id="single_cell_import",
            name="Single-Cell Data Import Wizard",
            description="Import single-cell RNA-seq data from various formats with automatic "
            "QC and preprocessing",
            input_formats=[".mtx", ".h5", ".h5ad", ".loom", ".csv", ".tsv"],
            output_formats=[".h5ad", ".rds", ".loom"],
        )

    def _initialize_parameters(self) -> List[ConversionParameter]:
        """Initialize scenario parameters"""
        return [
            # Input format
            ConversionParameter(
                name="input_format",
                label="Input Format",
                description="Single-cell data input format",
                type=ParameterType.SELECT,
                required=True,
                options=[
                    {
                        "value": SingleCellInputFormat.MTX_10X.value,
                        "label": "10X Genomics MTX (matrix.mtx + genes.tsv + barcodes.tsv)",
                    },
                    {
                        "value": SingleCellInputFormat.HDF5_10X.value,
                        "label": "10X Genomics HDF5 (.h5)",
                    },
                    {
                        "value": SingleCellInputFormat.H5AD.value,
                        "label": "AnnData (.h5ad)",
                    },
                    {
                        "value": SingleCellInputFormat.LOOM.value,
                        "label": "Loom (.loom)",
                    },
                    {
                        "value": SingleCellInputFormat.CSV.value,
                        "label": "CSV/TSV matrix",
                    },
                ],
            ),
            # Matrix directory (for 10X MTX)
            ConversionParameter(
                name="matrix_directory",
                label="Matrix Directory",
                description="Directory containing matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv",
                type=ParameterType.TEXT,
                required=False,
            ),
            # Species
            ConversionParameter(
                name="species",
                label="Species",
                description="Species of the sample",
                type=ParameterType.SELECT,
                required=True,
                default_value=Species.HUMAN.value,
                options=[
                    {"value": Species.HUMAN.value, "label": "Human (Homo sapiens)"},
                    {"value": Species.MOUSE.value, "label": "Mouse (Mus musculus)"},
                    {"value": Species.RAT.value, "label": "Rat (Rattus norvegicus)"},
                    {
                        "value": Species.ZEBRAFISH.value,
                        "label": "Zebrafish (Danio rerio)",
                    },
                    {
                        "value": Species.FLY.value,
                        "label": "Fruit fly (Drosophila melanogaster)",
                    },
                    {"value": Species.WORM.value, "label": "C. elegans"},
                ],
            ),
            # Genome version
            ConversionParameter(
                name="genome_version",
                label="Genome Version",
                description="Reference genome version",
                type=ParameterType.SELECT,
                required=True,
                default_value=GenomeVersion.GRCH38.value,
                options=[
                    {"value": GenomeVersion.GRCH38.value, "label": "GRCh38 (Human)"},
                    {"value": GenomeVersion.GRCH37.value, "label": "GRCh37 (Human)"},
                    {"value": GenomeVersion.GRCM39.value, "label": "GRCm39 (Mouse)"},
                    {"value": GenomeVersion.GRCM38.value, "label": "GRCm38 (Mouse)"},
                ],
            ),
            # QC: Minimum cells per gene
            ConversionParameter(
                name="min_cells",
                label="Minimum Cells per Gene",
                description="Filter genes detected in fewer than this many cells",
                type=ParameterType.NUMBER,
                required=True,
                default_value=3,
                validation_rules={"min": 0, "step": 1},
            ),
            # QC: Minimum genes per cell
            ConversionParameter(
                name="min_genes",
                label="Minimum Genes per Cell",
                description="Filter cells with fewer than this many detected genes",
                type=ParameterType.NUMBER,
                required=True,
                default_value=200,
                validation_rules={"min": 0, "step": 1},
            ),
            # QC: Maximum genes per cell
            ConversionParameter(
                name="max_genes",
                label="Maximum Genes per Cell",
                description="Filter cells with more than this many genes (doublets)",
                type=ParameterType.NUMBER,
                required=False,
                default_value=5000,
                validation_rules={"min": 0, "step": 1},
            ),
            # QC: Maximum mitochondrial percentage
            ConversionParameter(
                name="max_mito_percent",
                label="Maximum Mitochondrial %",
                description="Filter cells with more than this % of mitochondrial genes",
                type=ParameterType.NUMBER,
                required=True,
                default_value=10.0,
                validation_rules={"min": 0, "max": 100, "step": 0.1},
            ),
            # QC: Maximum ribosomal percentage (optional)
            ConversionParameter(
                name="max_ribo_percent",
                label="Maximum Ribosomal %",
                description="Filter cells with more than this % of ribosomal genes (optional)",
                type=ParameterType.NUMBER,
                required=False,
                validation_rules={"min": 0, "max": 100, "step": 0.1},
            ),
            # Normalization
            ConversionParameter(
                name="normalize",
                label="Normalize Data",
                description="Apply log-normalization (log1p after scaling to 10,000)",
                type=ParameterType.BOOLEAN,
                required=True,
                default_value=True,
            ),
            # Find variable features
            ConversionParameter(
                name="find_variable_features",
                label="Find Variable Features",
                description="Identify highly variable genes for downstream analysis",
                type=ParameterType.BOOLEAN,
                required=True,
                default_value=True,
            ),
            # Number of variable features
            ConversionParameter(
                name="n_variable_features",
                label="Number of Variable Features",
                description="Number of highly variable genes to identify",
                type=ParameterType.NUMBER,
                required=False,
                default_value=2000,
                validation_rules={"min": 100, "max": 10000, "step": 100},
            ),
            # Scale data
            ConversionParameter(
                name="scale_data",
                label="Scale Data",
                description="Scale expression values (z-score normalization)",
                type=ParameterType.BOOLEAN,
                required=True,
                default_value=False,
            ),
            # PCA
            ConversionParameter(
                name="run_pca",
                label="Run PCA",
                description="Run Principal Component Analysis",
                type=ParameterType.BOOLEAN,
                required=True,
                default_value=False,
            ),
            # Number of PCs
            ConversionParameter(
                name="n_pcs",
                label="Number of PCs",
                description="Number of principal components to compute",
                type=ParameterType.NUMBER,
                required=False,
                default_value=50,
                validation_rules={"min": 10, "max": 100, "step": 5},
            ),
            # Cell metadata file
            ConversionParameter(
                name="cell_metadata_file",
                label="Cell Metadata File",
                description="CSV/TSV file with cell annotations (optional)",
                type=ParameterType.FILE,
                required=False,
            ),
            # Sample name
            ConversionParameter(
                name="sample_name",
                label="Sample Name",
                description="Name for this sample (will be added to metadata)",
                type=ParameterType.TEXT,
                required=False,
                default_value="sample1",
            ),
            # Output format
            ConversionParameter(
                name="output_format",
                label="Output Format",
                description="Output file format",
                type=ParameterType.SELECT,
                required=True,
                default_value=OutputFormat.H5AD.value,
                options=[
                    {
                        "value": OutputFormat.H5AD.value,
                        "label": "AnnData h5ad (Python/scanpy)",
                    },
                    {"value": OutputFormat.RDS.value, "label": "Seurat RDS (R)"},
                    {"value": OutputFormat.LOOM.value, "label": "Loom format"},
                ],
            ),
        ]

    def detect_format(self, file_path: Path) -> float:
        """Detect if file is single-cell data"""
        try:
            suffix = file_path.suffix.lower()

            # Check file extension
            if suffix == ".h5ad":
                return 0.9
            elif suffix == ".loom":
                return 0.9
            elif suffix == ".h5":
                # Could be 10X HDF5
                return 0.7
            elif suffix == ".mtx":
                # Matrix Market format
                return 0.7
            elif suffix in [".csv", ".tsv", ".txt"]:
                # Check if looks like expression matrix
                delimiter = "\t" if suffix == ".tsv" else ","
                df = pd.read_csv(file_path, sep=delimiter, nrows=100)

                # Single-cell data typically has many cells (columns) and many genes (rows)
                if df.shape[1] > 100 and df.shape[0] > 1000:
                    return 0.5

                return 0.2

            return 0.0

        except Exception:
            return 0.0

    def validate_input(
        self, file_path: Path, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """Validate input file"""
        messages = []

        input_format = parameters.get("input_format")

        try:
            if input_format == SingleCellInputFormat.MTX_10X.value:
                # Check for required files
                matrix_dir = parameters.get("matrix_directory")
                if not matrix_dir:
                    messages.append(
                        ValidationMessage(
                            level=ValidationLevel.ERROR,
                            message="Matrix directory is required for 10X MTX format",
                        )
                    )
                else:
                    matrix_dir_path = Path(matrix_dir)
                    required_files = ["matrix.mtx", "barcodes.tsv", "genes.tsv"]

                    for req_file in required_files:
                        if not (matrix_dir_path / req_file).exists():
                            # Try alternative names
                            if (
                                req_file == "genes.tsv"
                                and (matrix_dir_path / "features.tsv").exists()
                            ):
                                continue
                            messages.append(
                                ValidationMessage(
                                    level=ValidationLevel.WARNING,
                                    message=f"Required file not found: {req_file}",
                                )
                            )

            elif input_format == SingleCellInputFormat.HDF5_10X.value:
                if file_path.suffix.lower() != ".h5":
                    messages.append(
                        ValidationMessage(
                            level=ValidationLevel.ERROR,
                            message="File must have .h5 extension for 10X HDF5 format",
                        )
                    )

            elif input_format == SingleCellInputFormat.H5AD.value:
                if file_path.suffix.lower() != ".h5ad":
                    messages.append(
                        ValidationMessage(
                            level=ValidationLevel.ERROR,
                            message="File must have .h5ad extension for AnnData format",
                        )
                    )

            elif input_format == SingleCellInputFormat.LOOM.value:
                if file_path.suffix.lower() != ".loom":
                    messages.append(
                        ValidationMessage(
                            level=ValidationLevel.ERROR,
                            message="File must have .loom extension for Loom format",
                        )
                    )

            # QC parameter validation
            min_cells = parameters.get("min_cells", 3)
            min_genes = parameters.get("min_genes", 200)
            max_genes = parameters.get("max_genes")

            if max_genes and max_genes < min_genes:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="Maximum genes must be greater than minimum genes",
                    )
                )

            messages.append(
                ValidationMessage(
                    level=ValidationLevel.INFO,
                    message=f"QC filters: min_cells={min_cells}, min_genes={min_genes}, "
                    f"max_mito={parameters.get('max_mito_percent', 10)}%",
                )
            )

        except Exception as e:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR,
                    message=f"Error validating file: {str(e)}",
                )
            )

        return messages

    def validate_parameters(
        self, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """Validate conversion parameters"""
        messages = []

        # Check QC thresholds
        min_genes = parameters.get("min_genes", 0)
        max_genes = parameters.get("max_genes")

        if max_genes and max_genes <= min_genes:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR,
                    message=f"max_genes ({max_genes}) must be greater than min_genes ({min_genes})",
                )
            )

        max_mito = parameters.get("max_mito_percent", 0)
        if max_mito < 0 or max_mito > 100:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR,
                    message=f"max_mito_percent must be between 0 and 100",
                )
            )

        # Check PCA parameters
        if parameters.get("run_pca"):
            n_pcs = parameters.get("n_pcs", 50)
            if n_pcs < 2:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR, message="n_pcs must be at least 2"
                    )
                )

        # Check variable features
        if parameters.get("find_variable_features"):
            n_var_features = parameters.get("n_variable_features", 2000)
            if n_var_features < 10:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.WARNING,
                        message=f"n_variable_features is very low ({n_var_features})",
                    )
                )

        return messages

    def generate_preview(
        self, file_path: Path, parameters: Dict[str, Any]
    ) -> ConversionPreview:
        """Generate conversion preview"""

        input_format = parameters.get("input_format")

        # Mock statistics (real implementation would load actual data)
        stats = {
            "input_format": input_format,
            "cells_before_qc": 5000,
            "genes_before_qc": 20000,
            "estimated_cells_after_qc": 4500,
            "estimated_genes_after_qc": 15000,
            "qc_filters": {
                "min_cells": parameters.get("min_cells", 3),
                "min_genes": parameters.get("min_genes", 200),
                "max_genes": parameters.get("max_genes"),
                "max_mito_percent": parameters.get("max_mito_percent", 10.0),
            },
            "preprocessing": {
                "normalize": parameters.get("normalize", True),
                "find_variable_features": parameters.get(
                    "find_variable_features", True
                ),
                "scale_data": parameters.get("scale_data", False),
                "run_pca": parameters.get("run_pca", False),
            },
        }

        messages = [
            ValidationMessage(
                level=ValidationLevel.INFO,
                message=f"Input: ~5,000 cells, ~20,000 genes",
            ),
            ValidationMessage(
                level=ValidationLevel.INFO,
                message=f"Expected after QC: ~4,500 cells, ~15,000 genes",
            ),
            ValidationMessage(
                level=ValidationLevel.INFO,
                message=f"Estimated filtering: ~10% cells, ~25% genes",
            ),
        ]

        # Preview sample data (mock)
        sample_data = {
            "columns": [
                "cell_barcode",
                "n_genes",
                "n_counts",
                "mito_percent",
                "pass_qc",
            ],
            "rows": [
                ["AAACCTGAGAAACCAT-1", 2500, 8000, 3.5, True],
                ["AAACCTGAGACAGACC-1", 1800, 5500, 5.2, True],
                ["AAACCTGAGGCCCGTT-1", 3200, 12000, 2.1, True],
                ["AAACCTGCAATCGAAA-1", 150, 400, 15.2, False],
                ["AAACCTGGTCCGTTAA-1", 2800, 9500, 4.8, True],
            ],
        }

        return ConversionPreview(
            sample_data=sample_data,
            validation_messages=messages,
            estimated_time_seconds=60,  # Rough estimate
            statistics=stats,
        )

    def convert(
        self,
        input_path: Path,
        output_path: Path,
        parameters: Dict[str, Any],
        progress_callback: Optional[Callable[[ConversionProgress], None]] = None,
    ) -> Dict[str, Any]:
        """Perform conversion"""

        def report_progress(stage: str, current: int, total: int):
            if progress_callback:
                progress_callback(
                    ConversionProgress(
                        stage=stage,
                        current_step=current,
                        total_steps=total,
                        message=f"{stage}...",
                    )
                )

        input_format = parameters.get("input_format")

        report_progress("Loading data", 1, 7)

        # Mock data loading (real implementation would use scanpy/Seurat)
        n_cells_raw = 5000
        n_genes_raw = 20000

        report_progress("Calculating QC metrics", 2, 7)

        # Mock QC metrics calculation

        report_progress("Filtering cells", 3, 7)

        # Apply QC filters
        min_genes = parameters.get("min_genes", 200)
        max_genes = parameters.get("max_genes", 5000)
        max_mito = parameters.get("max_mito_percent", 10.0)

        n_cells_after_qc = int(n_cells_raw * 0.9)  # Mock: 90% pass

        report_progress("Filtering genes", 4, 7)

        min_cells = parameters.get("min_cells", 3)
        n_genes_after_qc = int(n_genes_raw * 0.75)  # Mock: 75% pass

        report_progress("Normalization", 5, 7)

        # Apply normalization
        if parameters.get("normalize"):
            # log1p(counts / total_counts * 10000)
            pass

        report_progress("Finding variable features", 6, 7)

        if parameters.get("find_variable_features"):
            n_var_features = parameters.get("n_variable_features", 2000)
            # Identify highly variable genes
            pass

        report_progress("Saving output", 7, 7)

        output_format = parameters.get("output_format", OutputFormat.H5AD.value)

        # Save output (mock - would use scanpy/Seurat)
        # For now, just create a placeholder file
        output_path.write_text(
            json.dumps(
                {
                    "format": output_format,
                    "n_cells": n_cells_after_qc,
                    "n_genes": n_genes_after_qc,
                }
            )
        )

        return {
            "n_cells_raw": n_cells_raw,
            "n_genes_raw": n_genes_raw,
            "n_cells_after_qc": n_cells_after_qc,
            "n_genes_after_qc": n_genes_after_qc,
            "cells_filtered": n_cells_raw - n_cells_after_qc,
            "genes_filtered": n_genes_raw - n_genes_after_qc,
            "normalization_applied": parameters.get("normalize", True),
            "variable_features_found": parameters.get("find_variable_features", True),
            "output_format": output_format,
        }


# Singleton instance
_single_cell_import_scenario = None


def get_single_cell_import_scenario() -> SingleCellImportScenario:
    """Get singleton instance of single-cell import scenario"""
    global _single_cell_import_scenario
    if _single_cell_import_scenario is None:
        _single_cell_import_scenario = SingleCellImportScenario()
    return _single_cell_import_scenario
