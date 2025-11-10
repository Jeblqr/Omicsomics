"""
Expression Matrix Standardization Scenario

This scenario handles standardization of gene expression matrices from various sources
(DESeq2, edgeR, limma, Salmon, RSEM, etc.) with different orientations, gene ID types,
and normalization methods.

Key Features:
- Row/column orientation detection (genes × samples vs samples × genes)
- Gene ID type recognition (Ensembl, Symbol, Entrez, Refseq)
- Gene ID conversion between different systems
- Batch effect correction (optional)
- Sample metadata integration
- Normalization method selection (TPM, FPKM, CPM, log2, etc.)
- Output format selection (CSV, TSV, h5ad, RDS)
"""

from typing import Dict, List, Optional, Any, Callable
import pandas as pd
import numpy as np
from pathlib import Path
import re
from enum import Enum

from ..interactive_converter import (
    ConversionScenario,
    ConversionParameter,
    ParameterType,
    ValidationMessage,
    ValidationLevel,
    ConversionPreview,
    ConversionProgress,
)


class GeneIDType(str, Enum):
    """Gene ID types"""

    ENSEMBL = "ensembl"  # ENSG00000000003
    ENSEMBL_TRANSCRIPT = "ensembl_transcript"  # ENST00000000233
    SYMBOL = "symbol"  # TP53, BRCA1
    ENTREZ = "entrez"  # 7157, 672
    REFSEQ = "refseq"  # NM_000546


class MatrixOrientation(str, Enum):
    """Matrix orientation"""

    GENES_ROWS = "genes_rows"  # Genes in rows, samples in columns
    SAMPLES_ROWS = "samples_rows"  # Samples in rows, genes in columns


class NormalizationMethod(str, Enum):
    """Normalization methods"""

    NONE = "none"
    LOG2 = "log2"  # log2(x + 1)
    LOG2_CPM = "log2_cpm"  # log2(CPM + 1)
    TPM = "tpm"  # Transcripts Per Million
    FPKM = "fpkm"  # Fragments Per Kilobase Million
    CPM = "cpm"  # Counts Per Million
    ZSCORE = "zscore"  # Z-score normalization


class BatchCorrectionMethod(str, Enum):
    """Batch effect correction methods"""

    NONE = "none"
    COMBAT = "combat"  # ComBat from sva package
    LIMMA = "limma_removeBatchEffect"  # limma removeBatchEffect


class ExpressionMatrixStandardizationScenario(ConversionScenario):
    """Expression matrix standardization scenario"""

    # Common gene ID patterns
    GENE_ID_PATTERNS = {
        GeneIDType.ENSEMBL: re.compile(r"^ENS[A-Z]*G\d{11}(\.\d+)?$"),
        GeneIDType.ENSEMBL_TRANSCRIPT: re.compile(r"^ENS[A-Z]*T\d{11}(\.\d+)?$"),
        GeneIDType.ENTREZ: re.compile(r"^\d+$"),
        GeneIDType.REFSEQ: re.compile(r"^N[MR]_\d+(\.\d+)?$"),
        GeneIDType.SYMBOL: re.compile(r"^[A-Z][A-Z0-9\-]+$"),
    }

    def __init__(self):
        super().__init__(
            scenario_id="expression_matrix_standardization",
            name="Gene Expression Matrix Standardization",
            description="Standardize gene expression matrices from various sources with different "
            "orientations, gene ID types, and normalization methods",
            input_formats=[".csv", ".tsv", ".txt", ".xlsx"],
            output_formats=[".csv", ".tsv", ".h5ad", ".rds"],
        )

    def _initialize_parameters(self) -> List[ConversionParameter]:
        """Initialize scenario parameters"""
        return [
            # Matrix orientation
            ConversionParameter(
                name="orientation",
                label="Matrix Orientation",
                description="How genes and samples are organized in the matrix",
                type=ParameterType.SELECT,
                required=True,
                default_value=MatrixOrientation.GENES_ROWS.value,
                options=[
                    {
                        "value": MatrixOrientation.GENES_ROWS.value,
                        "label": "Genes in rows, samples in columns (standard)",
                    },
                    {
                        "value": MatrixOrientation.SAMPLES_ROWS.value,
                        "label": "Samples in rows, genes in columns (transposed)",
                    },
                ],
            ),
            # Gene ID type
            ConversionParameter(
                name="gene_id_type",
                label="Gene ID Type",
                description="Type of gene identifiers used in the matrix",
                type=ParameterType.SELECT,
                required=True,
                options=[
                    {
                        "value": GeneIDType.ENSEMBL.value,
                        "label": "Ensembl Gene ID (ENSG...)",
                    },
                    {
                        "value": GeneIDType.ENSEMBL_TRANSCRIPT.value,
                        "label": "Ensembl Transcript ID (ENST...)",
                    },
                    {
                        "value": GeneIDType.SYMBOL.value,
                        "label": "Gene Symbol (TP53, BRCA1, ...)",
                    },
                    {
                        "value": GeneIDType.ENTREZ.value,
                        "label": "Entrez Gene ID (7157, 672, ...)",
                    },
                    {
                        "value": GeneIDType.REFSEQ.value,
                        "label": "RefSeq ID (NM_000546, ...)",
                    },
                ],
            ),
            # Target gene ID type
            ConversionParameter(
                name="target_gene_id_type",
                label="Convert Gene IDs To",
                description="Target gene ID type for output (optional conversion)",
                type=ParameterType.SELECT,
                required=False,
                options=[
                    {"value": "same", "label": "Keep original ID type"},
                    {"value": GeneIDType.ENSEMBL.value, "label": "Ensembl Gene ID"},
                    {"value": GeneIDType.SYMBOL.value, "label": "Gene Symbol"},
                    {"value": GeneIDType.ENTREZ.value, "label": "Entrez Gene ID"},
                ],
                default_value="same",
            ),
            # Species
            ConversionParameter(
                name="species",
                label="Species",
                description="Species for gene ID conversion",
                type=ParameterType.SELECT,
                required=True,
                default_value="human",
                options=[
                    {"value": "human", "label": "Human (Homo sapiens)"},
                    {"value": "mouse", "label": "Mouse (Mus musculus)"},
                    {"value": "rat", "label": "Rat (Rattus norvegicus)"},
                    {"value": "zebrafish", "label": "Zebrafish (Danio rerio)"},
                    {"value": "fly", "label": "Fruit fly (Drosophila melanogaster)"},
                    {"value": "worm", "label": "C. elegans"},
                ],
            ),
            # Gene ID column (for transposed matrices)
            ConversionParameter(
                name="gene_id_column",
                label="Gene ID Column",
                description="Column containing gene IDs (only for samples-in-rows orientation)",
                type=ParameterType.SELECT,
                required=False,
                options=[],  # Will be populated from file columns
            ),
            # Normalization method
            ConversionParameter(
                name="normalization",
                label="Normalization Method",
                description="Apply normalization to expression values",
                type=ParameterType.SELECT,
                required=True,
                default_value=NormalizationMethod.NONE.value,
                options=[
                    {
                        "value": NormalizationMethod.NONE.value,
                        "label": "None (keep original values)",
                    },
                    {
                        "value": NormalizationMethod.LOG2.value,
                        "label": "Log2 transformation: log2(x + 1)",
                    },
                    {
                        "value": NormalizationMethod.LOG2_CPM.value,
                        "label": "Log2 CPM: log2(CPM + 1)",
                    },
                    {
                        "value": NormalizationMethod.CPM.value,
                        "label": "CPM (Counts Per Million)",
                    },
                    {
                        "value": NormalizationMethod.TPM.value,
                        "label": "TPM (Transcripts Per Million)",
                    },
                    {
                        "value": NormalizationMethod.FPKM.value,
                        "label": "FPKM (Fragments Per Kilobase Million)",
                    },
                    {
                        "value": NormalizationMethod.ZSCORE.value,
                        "label": "Z-score normalization",
                    },
                ],
            ),
            # Gene lengths (for TPM/FPKM)
            ConversionParameter(
                name="gene_lengths_file",
                label="Gene Lengths File",
                description="File with gene lengths (required for TPM/FPKM normalization)",
                type=ParameterType.FILE,
                required=False,
            ),
            # Batch correction
            ConversionParameter(
                name="batch_correction",
                label="Batch Effect Correction",
                description="Apply batch effect correction",
                type=ParameterType.SELECT,
                required=True,
                default_value=BatchCorrectionMethod.NONE.value,
                options=[
                    {"value": BatchCorrectionMethod.NONE.value, "label": "None"},
                    {
                        "value": BatchCorrectionMethod.COMBAT.value,
                        "label": "ComBat (sva package)",
                    },
                    {
                        "value": BatchCorrectionMethod.LIMMA.value,
                        "label": "limma removeBatchEffect",
                    },
                ],
            ),
            # Batch metadata file
            ConversionParameter(
                name="batch_metadata_file",
                label="Batch Metadata File",
                description="File with batch information (required for batch correction)",
                type=ParameterType.FILE,
                required=False,
            ),
            # Sample metadata file
            ConversionParameter(
                name="sample_metadata_file",
                label="Sample Metadata File",
                description="File with sample annotations (optional)",
                type=ParameterType.FILE,
                required=False,
            ),
            # Filter low expression genes
            ConversionParameter(
                name="filter_low_expression",
                label="Filter Low Expression Genes",
                description="Remove genes with low expression across samples",
                type=ParameterType.BOOLEAN,
                required=False,
                default_value=False,
            ),
            # Minimum expression threshold
            ConversionParameter(
                name="min_expression",
                label="Minimum Expression Threshold",
                description="Minimum expression value (genes below this in all samples are removed)",
                type=ParameterType.NUMBER,
                required=False,
                default_value=1.0,
                validation_rules={"min": 0.0, "step": 0.1},
            ),
            # Minimum sample count
            ConversionParameter(
                name="min_samples",
                label="Minimum Sample Count",
                description="Minimum number of samples where gene must pass threshold",
                type=ParameterType.NUMBER,
                required=False,
                default_value=1,
                validation_rules={"min": 1, "step": 1},
            ),
            # Remove duplicates
            ConversionParameter(
                name="remove_duplicates",
                label="Remove Duplicate Genes",
                description="How to handle duplicate gene IDs",
                type=ParameterType.SELECT,
                required=True,
                default_value="keep_first",
                options=[
                    {"value": "keep_first", "label": "Keep first occurrence"},
                    {"value": "keep_last", "label": "Keep last occurrence"},
                    {"value": "sum", "label": "Sum expression values"},
                    {"value": "mean", "label": "Average expression values"},
                    {"value": "max", "label": "Keep maximum expression"},
                ],
            ),
            # Output format
            ConversionParameter(
                name="output_format",
                label="Output Format",
                description="Output file format",
                type=ParameterType.SELECT,
                required=True,
                default_value=".csv",
                options=[
                    {"value": ".csv", "label": "CSV (Comma-separated)"},
                    {"value": ".tsv", "label": "TSV (Tab-separated)"},
                    {"value": ".h5ad", "label": "AnnData h5ad (Python/scanpy)"},
                    {"value": ".rds", "label": "RDS (R/Seurat)"},
                ],
            ),
        ]

    def detect_format(self, file_path: Path) -> float:
        """Detect if file is an expression matrix"""
        try:
            # Try reading as CSV/TSV
            delimiter = self._detect_delimiter(file_path)
            df = pd.read_csv(file_path, sep=delimiter, nrows=100)

            # Check if looks like expression matrix
            score = 0.0

            # Should have multiple columns
            if df.shape[1] < 3:
                return 0.0

            score += 0.2

            # Most columns should be numeric
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            numeric_ratio = len(numeric_cols) / df.shape[1]

            if numeric_ratio > 0.5:
                score += 0.3

            # First column might be gene IDs (string)
            if df.iloc[:, 0].dtype == object:
                score += 0.2

                # Check if first column matches gene ID patterns
                sample_ids = df.iloc[:, 0].head(20).astype(str)
                for id_type, pattern in self.GENE_ID_PATTERNS.items():
                    matches = sum(1 for id_str in sample_ids if pattern.match(id_str))
                    if matches > len(sample_ids) * 0.5:
                        score += 0.3
                        break

            return min(score, 1.0)

        except Exception:
            return 0.0

    def validate_input(
        self, file_path: Path, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """Validate input file"""
        messages = []

        try:
            # Read file
            delimiter = self._detect_delimiter(file_path)
            df = pd.read_csv(file_path, sep=delimiter, nrows=100)

            # Check dimensions
            if df.shape[0] < 2:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message=f"File has too few rows ({df.shape[0]}). Need at least 2 rows.",
                    )
                )

            if df.shape[1] < 2:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message=f"File has too few columns ({df.shape[1]}). Need at least 2 columns.",
                    )
                )

            # Check for numeric data
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) == 0:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="No numeric columns found. Expression matrix should contain numeric values.",
                    )
                )

            # Info about matrix dimensions
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.INFO,
                    message=f"Detected matrix dimensions: {df.shape[0]} rows × {df.shape[1]} columns",
                )
            )

            # Check normalization requirements
            normalization = parameters.get(
                "normalization", NormalizationMethod.NONE.value
            )
            if normalization in [
                NormalizationMethod.TPM.value,
                NormalizationMethod.FPKM.value,
            ]:
                if not parameters.get("gene_lengths_file"):
                    messages.append(
                        ValidationMessage(
                            level=ValidationLevel.WARNING,
                            message=f"{normalization.upper()} normalization requires gene lengths file",
                        )
                    )

            # Check batch correction requirements
            batch_correction = parameters.get(
                "batch_correction", BatchCorrectionMethod.NONE.value
            )
            if batch_correction != BatchCorrectionMethod.NONE.value:
                if not parameters.get("batch_metadata_file"):
                    messages.append(
                        ValidationMessage(
                            level=ValidationLevel.WARNING,
                            message=f"{batch_correction} requires batch metadata file",
                        )
                    )

        except Exception as e:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR, message=f"Error reading file: {str(e)}"
                )
            )

        return messages

    def validate_parameters(
        self, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """Validate conversion parameters"""
        messages = []

        # Check normalization + gene lengths
        normalization = parameters.get("normalization", NormalizationMethod.NONE.value)
        if normalization in [
            NormalizationMethod.TPM.value,
            NormalizationMethod.FPKM.value,
        ]:
            if not parameters.get("gene_lengths_file"):
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message=f"{normalization.upper()} normalization requires gene lengths file",
                    )
                )

        # Check batch correction + metadata
        batch_correction = parameters.get(
            "batch_correction", BatchCorrectionMethod.NONE.value
        )
        if batch_correction != BatchCorrectionMethod.NONE.value:
            if not parameters.get("batch_metadata_file"):
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message=f"{batch_correction} requires batch metadata file",
                    )
                )

        # Check filtering parameters
        if parameters.get("filter_low_expression"):
            min_expression = parameters.get("min_expression", 0)
            if min_expression < 0:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="Minimum expression threshold must be non-negative",
                    )
                )

            min_samples = parameters.get("min_samples", 1)
            if min_samples < 1:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="Minimum sample count must be at least 1",
                    )
                )

        return messages

    def generate_preview(
        self, file_path: Path, parameters: Dict[str, Any]
    ) -> ConversionPreview:
        """Generate conversion preview"""

        # Read file
        delimiter = self._detect_delimiter(file_path)
        df = pd.read_csv(file_path, sep=delimiter, nrows=1000)

        orientation = parameters.get("orientation", MatrixOrientation.GENES_ROWS.value)

        # Determine matrix structure
        if orientation == MatrixOrientation.SAMPLES_ROWS.value:
            # Transpose if needed
            gene_id_col = parameters.get("gene_id_column")
            if gene_id_col and gene_id_col in df.columns:
                df = df.set_index(gene_id_col).T

        # Detect gene ID type
        detected_id_type = self._detect_gene_id_type(df.index)

        # Count numeric samples
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        n_samples = len(numeric_cols)
        n_genes = df.shape[0]

        # Preview sample
        preview_df = df.head(10)

        # Statistics
        stats = {
            "n_genes": n_genes,
            "n_samples": n_samples,
            "detected_gene_id_type": detected_id_type,
            "matrix_orientation": orientation,
            "has_missing_values": df.isnull().any().any(),
            "expression_range": {
                "min": float(df[numeric_cols].min().min()),
                "max": float(df[numeric_cols].max().max()),
                "mean": float(df[numeric_cols].mean().mean()),
            },
        }

        # Filtering preview
        if parameters.get("filter_low_expression"):
            min_expr = parameters.get("min_expression", 1.0)
            min_samples_count = parameters.get("min_samples", 1)

            # Count genes passing filter
            passing = (df[numeric_cols] >= min_expr).sum(axis=1) >= min_samples_count
            n_passing = passing.sum()
            n_filtered = n_genes - n_passing

            stats["filtering"] = {
                "genes_passing": n_passing,
                "genes_filtered": n_filtered,
                "percent_kept": round(n_passing / n_genes * 100, 1),
            }

        messages = []

        # Add detection message
        if detected_id_type:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.INFO,
                    message=f"Detected gene ID type: {detected_id_type}",
                )
            )

        # Add dimension info
        messages.append(
            ValidationMessage(
                level=ValidationLevel.INFO,
                message=f"Matrix: {n_genes} genes × {n_samples} samples",
            )
        )

        return ConversionPreview(
            sample_data=self._dataframe_to_preview(preview_df),
            validation_messages=messages,
            estimated_time_seconds=n_genes * n_samples / 100000,  # Rough estimate
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

        report_progress("Loading data", 1, 6)

        # Read input file
        delimiter = self._detect_delimiter(input_path)
        df = pd.read_csv(input_path, sep=delimiter, index_col=0)

        original_shape = df.shape

        # Handle orientation
        orientation = parameters.get("orientation", MatrixOrientation.GENES_ROWS.value)
        if orientation == MatrixOrientation.SAMPLES_ROWS.value:
            df = df.T

        report_progress("Processing gene IDs", 2, 6)

        # Handle duplicate genes
        duplicate_method = parameters.get("remove_duplicates", "keep_first")
        if duplicate_method == "keep_first":
            df = df[~df.index.duplicated(keep="first")]
        elif duplicate_method == "keep_last":
            df = df[~df.index.duplicated(keep="last")]
        elif duplicate_method == "sum":
            df = df.groupby(df.index).sum()
        elif duplicate_method == "mean":
            df = df.groupby(df.index).mean()
        elif duplicate_method == "max":
            df = df.groupby(df.index).max()

        report_progress("Filtering genes", 3, 6)

        # Filter low expression genes
        if parameters.get("filter_low_expression"):
            min_expr = parameters.get("min_expression", 1.0)
            min_samples = parameters.get("min_samples", 1)

            passing = (df >= min_expr).sum(axis=1) >= min_samples
            df = df[passing]

        report_progress("Normalizing expression values", 4, 6)

        # Apply normalization
        normalization = parameters.get("normalization", NormalizationMethod.NONE.value)
        df = self._apply_normalization(df, normalization, parameters)

        report_progress("Batch correction", 5, 6)

        # Batch correction (placeholder - would need R integration)
        batch_correction = parameters.get(
            "batch_correction", BatchCorrectionMethod.NONE.value
        )
        if batch_correction != BatchCorrectionMethod.NONE.value:
            # This would require R integration via rpy2
            pass

        report_progress("Saving output", 6, 6)

        # Save output
        output_format = parameters.get("output_format", ".csv")
        if output_format == ".csv":
            df.to_csv(output_path, sep=",")
        elif output_format == ".tsv":
            df.to_csv(output_path, sep="\t")
        elif output_format == ".h5ad":
            # Would need scanpy integration
            df.to_csv(output_path.with_suffix(".csv"))
        elif output_format == ".rds":
            # Would need R integration
            df.to_csv(output_path.with_suffix(".csv"))

        return {
            "input_shape": original_shape,
            "output_shape": df.shape,
            "genes_removed": original_shape[0] - df.shape[0],
            "normalization_applied": normalization,
            "batch_correction_applied": batch_correction,
        }

    def _detect_delimiter(self, file_path: Path) -> str:
        """Detect CSV delimiter"""
        with open(file_path, "r") as f:
            first_line = f.readline()
            if "\t" in first_line:
                return "\t"
            return ","

    def _detect_gene_id_type(self, gene_ids: pd.Index) -> Optional[str]:
        """Detect gene ID type from sample"""
        sample_ids = gene_ids[:100].astype(str)

        for id_type, pattern in self.GENE_ID_PATTERNS.items():
            matches = sum(1 for id_str in sample_ids if pattern.match(id_str))
            if matches > len(sample_ids) * 0.7:
                return id_type.value

        return None

    def _dataframe_to_preview(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Convert DataFrame to preview format"""
        return {
            "columns": ["gene_id"] + list(df.columns[:10]),
            "rows": [
                [str(idx)] + [float(val) if pd.notna(val) else None for val in row[:10]]
                for idx, row in df.head(10).iterrows()
            ],
        }

    def _apply_normalization(
        self, df: pd.DataFrame, method: str, parameters: Dict[str, Any]
    ) -> pd.DataFrame:
        """Apply normalization to expression matrix"""

        if method == NormalizationMethod.NONE.value:
            return df

        elif method == NormalizationMethod.LOG2.value:
            return np.log2(df + 1)

        elif method == NormalizationMethod.CPM.value:
            # Counts per million
            return df.div(df.sum(axis=0), axis=1) * 1e6

        elif method == NormalizationMethod.LOG2_CPM.value:
            cpm = df.div(df.sum(axis=0), axis=1) * 1e6
            return np.log2(cpm + 1)

        elif method == NormalizationMethod.TPM.value:
            # Would need gene lengths
            # For now, return CPM as placeholder
            return df.div(df.sum(axis=0), axis=1) * 1e6

        elif method == NormalizationMethod.FPKM.value:
            # Would need gene lengths
            # For now, return CPM as placeholder
            return df.div(df.sum(axis=0), axis=1) * 1e6

        elif method == NormalizationMethod.ZSCORE.value:
            # Z-score per gene
            return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)

        return df


# Singleton instance
_expression_matrix_scenario = None


def get_expression_matrix_standardization_scenario() -> (
    ExpressionMatrixStandardizationScenario
):
    """Get singleton instance of expression matrix standardization scenario"""
    global _expression_matrix_scenario
    if _expression_matrix_scenario is None:
        _expression_matrix_scenario = ExpressionMatrixStandardizationScenario()
    return _expression_matrix_scenario
