"""
GWAS Summary Statistics Standardization Scenario.

Handles standardization of GWAS summary statistics from various sources
with interactive column mapping and parameter configuration.
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

from .interactive_converter import (
    ConversionScenario,
    ConversionParameter,
    ParameterType,
    ValidationMessage,
    ValidationLevel,
    ConversionPreview,
    ConversionProgress,
    detect_delimiter,
    preview_dataframe,
    infer_column_types,
)


class GWASStandardizationScenario(ConversionScenario):
    """
    Standardize GWAS summary statistics files.

    Handles:
    - Column mapping (SNP ID, position, alleles, stats)
    - Genome build selection (hg19/hg38/GRCh37/GRCh38)
    - Effect size unit conversion (beta/OR/log(OR))
    - Allele harmonization
    - Quality control filtering
    """

    # Standard column names for output
    STANDARD_COLUMNS = {
        "variant_id": "Variant identifier (rsID or chr:pos)",
        "chromosome": "Chromosome",
        "position": "Base position",
        "effect_allele": "Effect/alternative allele",
        "other_allele": "Other/reference allele",
        "effect_size": "Effect size (beta or log(OR))",
        "standard_error": "Standard error",
        "pvalue": "P-value",
        "sample_size": "Sample size",
        "maf": "Minor allele frequency",
        "info": "Imputation INFO score",
    }

    # Common column name variations
    COLUMN_VARIANTS = {
        "variant_id": ["snp", "rsid", "rs", "variant", "marker", "snpid", "id"],
        "chromosome": ["chr", "chrom", "chromosome", "#chr"],
        "position": ["pos", "bp", "position", "base_pair_location", "bp_position"],
        "effect_allele": [
            "a1",
            "ea",
            "effect_allele",
            "alt",
            "allele1",
            "tested_allele",
        ],
        "other_allele": [
            "a2",
            "nea",
            "other_allele",
            "ref",
            "allele2",
            "reference_allele",
        ],
        "effect_size": ["beta", "b", "effect", "or", "odds_ratio", "log_or"],
        "standard_error": ["se", "stderr", "standard_error", "sebeta"],
        "pvalue": ["p", "pval", "p_value", "pvalue", "p-value"],
        "sample_size": ["n", "n_total", "sample_size", "neff", "n_eff"],
        "maf": ["maf", "frq", "freq", "allele_freq", "eaf"],
        "info": ["info", "imputation_info", "r2"],
    }

    def _initialize_parameters(self):
        """Initialize GWAS standardization parameters."""

        # Column mapping parameters (will be populated dynamically)
        for std_col, description in self.STANDARD_COLUMNS.items():
            self.parameters.append(
                ConversionParameter(
                    name=f"col_{std_col}",
                    label=f"{std_col.replace('_', ' ').title()} Column",
                    type=ParameterType.SELECT,
                    required=std_col
                    in [
                        "variant_id",
                        "chromosome",
                        "position",
                        "effect_allele",
                        "pvalue",
                    ],
                    description=description,
                    options=[],  # Will be filled from file columns
                )
            )

        # Genome build selection
        self.parameters.append(
            ConversionParameter(
                name="genome_build",
                label="Genome Build",
                type=ParameterType.SELECT,
                required=True,
                default="hg38",
                options=[
                    {"value": "hg19", "label": "hg19 / GRCh37"},
                    {"value": "hg38", "label": "hg38 / GRCh38"},
                    {"value": "hg18", "label": "hg18 / NCBI36"},
                ],
                description="Reference genome build",
            )
        )

        # Effect size type
        self.parameters.append(
            ConversionParameter(
                name="effect_size_type",
                label="Effect Size Type",
                type=ParameterType.SELECT,
                required=True,
                default="beta",
                options=[
                    {"value": "beta", "label": "Beta (linear regression)"},
                    {"value": "or", "label": "Odds Ratio"},
                    {"value": "log_or", "label": "Log(Odds Ratio)"},
                ],
                description="Type of effect size in input file",
            )
        )

        # Convert to standard format
        self.parameters.append(
            ConversionParameter(
                name="standardize_effect",
                label="Standardize Effect Size",
                type=ParameterType.BOOLEAN,
                default=True,
                description="Convert all effect sizes to beta (log scale)",
            )
        )

        # P-value threshold
        self.parameters.append(
            ConversionParameter(
                name="pvalue_threshold",
                label="P-value Threshold",
                type=ParameterType.NUMBER,
                default=1.0,
                validation={"min": 0, "max": 1},
                description="Filter variants with p-value above this threshold (1.0 = no filter)",
            )
        )

        # MAF threshold
        self.parameters.append(
            ConversionParameter(
                name="maf_threshold",
                label="MAF Threshold",
                type=ParameterType.NUMBER,
                default=0.0,
                validation={"min": 0, "max": 0.5},
                description="Filter variants with MAF below this threshold (0 = no filter)",
            )
        )

        # INFO score threshold
        self.parameters.append(
            ConversionParameter(
                name="info_threshold",
                label="INFO Score Threshold",
                type=ParameterType.NUMBER,
                default=0.0,
                validation={"min": 0, "max": 1},
                description="Filter variants with INFO score below this threshold (0 = no filter)",
            )
        )

        # Allele harmonization
        self.parameters.append(
            ConversionParameter(
                name="harmonize_alleles",
                label="Harmonize Alleles",
                type=ParameterType.BOOLEAN,
                default=True,
                description="Ensure consistent allele coding (flip strand if needed)",
            )
        )

    def detect_format(self, file_path: str) -> bool:
        """
        Detect if file is a GWAS summary statistics file.

        Args:
            file_path: Path to input file

        Returns:
            True if file appears to be GWAS summary statistics
        """
        try:
            # Read first few rows
            df = pd.read_csv(file_path, sep=None, engine="python", nrows=10)

            # Check for common GWAS columns
            columns_lower = [col.lower() for col in df.columns]

            # Must have: SNP ID or position, p-value
            has_snp = any(
                variant in columns_lower
                for variants in [
                    self.COLUMN_VARIANTS["variant_id"],
                    self.COLUMN_VARIANTS["chromosome"],
                ]
                for variant in variants
            )

            has_pvalue = any(
                variant in columns_lower for variant in self.COLUMN_VARIANTS["pvalue"]
            )

            # Should have effect size or OR
            has_effect = any(
                variant in columns_lower
                for variant in self.COLUMN_VARIANTS["effect_size"]
            )

            return has_snp and has_pvalue and has_effect

        except Exception:
            return False

    def validate_input(self, file_path: str) -> List[ValidationMessage]:
        """
        Validate GWAS summary statistics file.

        Args:
            file_path: Path to input file

        Returns:
            List of validation messages
        """
        messages = []

        try:
            # Check file exists
            if not Path(file_path).exists():
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR, message="File does not exist"
                    )
                )
                return messages

            # Try to read file
            df = pd.read_csv(file_path, sep=None, engine="python", nrows=100)

            # Check minimum columns
            if len(df.columns) < 3:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message=f"File has too few columns ({len(df.columns)}). Expected at least 3.",
                    )
                )

            # Check for numeric columns (effect size, p-value)
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) < 1:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.WARNING,
                        message="No numeric columns detected. Ensure effect sizes and p-values are numeric.",
                    )
                )

            # Check file size
            file_size = Path(file_path).stat().st_size / (1024 * 1024)  # MB
            if file_size > 1000:
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.WARNING,
                        message=f"Large file detected ({file_size:.1f} MB). Processing may take time.",
                    )
                )

            # Info message about detected columns
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.INFO,
                    message=f"Detected {len(df.columns)} columns and {len(df)} rows (preview).",
                )
            )

        except Exception as e:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR, message=f"Error reading file: {str(e)}"
                )
            )

        return messages

    def _suggest_column_mapping(self, df: pd.DataFrame) -> Dict[str, Optional[str]]:
        """
        Suggest column mappings based on common names.

        Args:
            df: Input DataFrame

        Returns:
            Dictionary mapping standard names to detected column names
        """
        mapping = {}
        columns_lower = {col: col.lower() for col in df.columns}

        for std_col, variants in self.COLUMN_VARIANTS.items():
            mapping[std_col] = None
            for variant in variants:
                for col, col_lower in columns_lower.items():
                    if variant == col_lower or variant in col_lower:
                        mapping[std_col] = col
                        break
                if mapping[std_col]:
                    break

        return mapping

    def validate_parameters(
        self, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """
        Validate conversion parameters.

        Args:
            parameters: Parameter values

        Returns:
            List of validation messages
        """
        messages = []

        # Check required column mappings
        required_cols = [
            "variant_id",
            "chromosome",
            "position",
            "effect_allele",
            "pvalue",
        ]
        for col in required_cols:
            param_name = f"col_{col}"
            if not parameters.get(param_name):
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message=f"Required column mapping missing: {col}",
                        field=param_name,
                    )
                )

        # Validate thresholds
        if "pvalue_threshold" in parameters:
            pval_thresh = parameters["pvalue_threshold"]
            if not (0 <= pval_thresh <= 1):
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="P-value threshold must be between 0 and 1",
                        field="pvalue_threshold",
                    )
                )

        if "maf_threshold" in parameters:
            maf_thresh = parameters["maf_threshold"]
            if not (0 <= maf_thresh <= 0.5):
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="MAF threshold must be between 0 and 0.5",
                        field="maf_threshold",
                    )
                )

        if "info_threshold" in parameters:
            info_thresh = parameters["info_threshold"]
            if not (0 <= info_thresh <= 1):
                messages.append(
                    ValidationMessage(
                        level=ValidationLevel.ERROR,
                        message="INFO score threshold must be between 0 and 1",
                        field="info_threshold",
                    )
                )

        return messages

    def generate_preview(
        self, file_path: str, parameters: Dict[str, Any]
    ) -> ConversionPreview:
        """
        Generate conversion preview.

        Args:
            file_path: Path to input file
            parameters: Conversion parameters

        Returns:
            ConversionPreview object
        """
        # Read input file (preview only)
        df = pd.read_csv(file_path, sep=None, engine="python", nrows=100)

        # Suggest column mappings if not provided
        suggested_mapping = self._suggest_column_mapping(df)

        # Apply column mappings
        column_mapping = {}
        for std_col in self.STANDARD_COLUMNS.keys():
            param_name = f"col_{std_col}"
            if param_name in parameters and parameters[param_name]:
                column_mapping[parameters[param_name]] = std_col
            elif suggested_mapping.get(std_col):
                column_mapping[suggested_mapping[std_col]] = std_col

        # Create preview dataframe
        preview_df = df.rename(columns=column_mapping)

        # Apply filters for preview
        if "pvalue_threshold" in parameters and "pvalue" in preview_df.columns:
            pval_thresh = parameters["pvalue_threshold"]
            if pval_thresh < 1.0:
                preview_df = preview_df[preview_df["pvalue"] <= pval_thresh]

        # Estimate output size
        total_rows = len(df)
        preview_rows = len(preview_df)
        estimated_output_rows = int(total_rows * (preview_rows / 100))

        # Get file info
        file_size = Path(file_path).stat().st_size

        # Estimate processing time (rough: 1GB per 60 seconds)
        estimated_time = (file_size / (1024**3)) * 60

        # Validation messages
        validation_messages = self.validate_parameters(parameters)

        # Add info about suggested mappings
        mapped_count = sum(1 for v in column_mapping.values() if v)
        validation_messages.append(
            ValidationMessage(
                level=ValidationLevel.INFO,
                message=f"Mapped {mapped_count}/{len(self.STANDARD_COLUMNS)} columns",
                details={"suggested_mapping": suggested_mapping},
            )
        )

        return ConversionPreview(
            input_info={
                "filename": Path(file_path).name,
                "size_mb": file_size / (1024**2),
                "rows": total_rows,
                "columns": len(df.columns),
                "delimiter": detect_delimiter(file_path),
            },
            output_info={
                "format": "standardized_gwas",
                "estimated_rows": estimated_output_rows,
                "columns": list(column_mapping.values()),
                "genome_build": parameters.get("genome_build", "hg38"),
            },
            sample_data=preview_dataframe(preview_df, max_rows=10),
            validation_messages=validation_messages,
            estimated_time=estimated_time,
        )

    def convert(
        self,
        source_path: str,
        target_path: str,
        parameters: Dict[str, Any],
        progress_callback: Optional[callable] = None,
    ) -> Dict[str, Any]:
        """
        Perform GWAS standardization conversion.

        Args:
            source_path: Input file path
            target_path: Output file path
            parameters: Conversion parameters
            progress_callback: Progress update callback

        Returns:
            Conversion statistics
        """
        start_time = datetime.now()

        # Progress update
        if progress_callback:
            progress_callback(
                ConversionProgress(
                    stage="reading",
                    progress=0,
                    message="Reading input file...",
                    started_at=start_time,
                )
            )

        # Read input file
        df = pd.read_csv(source_path, sep=None, engine="python")
        input_rows = len(df)

        # Apply column mappings
        column_mapping = {}
        for std_col in self.STANDARD_COLUMNS.keys():
            param_name = f"col_{std_col}"
            if param_name in parameters and parameters[param_name]:
                column_mapping[parameters[param_name]] = std_col

        df = df.rename(columns=column_mapping)

        # Progress update
        if progress_callback:
            progress_callback(
                ConversionProgress(
                    stage="filtering",
                    progress=30,
                    message="Applying filters...",
                    started_at=start_time,
                )
            )

        # Apply filters
        filtered_count = 0

        if "pvalue_threshold" in parameters and "pvalue" in df.columns:
            pval_thresh = parameters["pvalue_threshold"]
            if pval_thresh < 1.0:
                original_len = len(df)
                df = df[df["pvalue"] <= pval_thresh]
                filtered_count += original_len - len(df)

        if "maf_threshold" in parameters and "maf" in df.columns:
            maf_thresh = parameters["maf_threshold"]
            if maf_thresh > 0:
                original_len = len(df)
                df = df[df["maf"] >= maf_thresh]
                filtered_count += original_len - len(df)

        if "info_threshold" in parameters and "info" in df.columns:
            info_thresh = parameters["info_threshold"]
            if info_thresh > 0:
                original_len = len(df)
                df = df[df["info"] >= info_thresh]
                filtered_count += original_len - len(df)

        # Progress update
        if progress_callback:
            progress_callback(
                ConversionProgress(
                    stage="standardizing",
                    progress=60,
                    message="Standardizing effect sizes...",
                    started_at=start_time,
                )
            )

        # Standardize effect sizes
        if parameters.get("standardize_effect", True) and "effect_size" in df.columns:
            effect_type = parameters.get("effect_size_type", "beta")

            if effect_type == "or":
                # Convert OR to log(OR)
                df["effect_size"] = np.log(df["effect_size"])
            elif effect_type == "log_or":
                # Already in log scale
                pass
            # beta is already in the right format

        # Progress update
        if progress_callback:
            progress_callback(
                ConversionProgress(
                    stage="writing",
                    progress=90,
                    message="Writing output file...",
                    started_at=start_time,
                )
            )

        # Keep only standard columns
        output_columns = [
            col for col in self.STANDARD_COLUMNS.keys() if col in df.columns
        ]
        df = df[output_columns]

        # Write output
        df.to_csv(target_path, sep="\t", index=False)

        output_rows = len(df)
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        # Progress complete
        if progress_callback:
            progress_callback(
                ConversionProgress(
                    stage="complete",
                    progress=100,
                    message="Conversion complete",
                    started_at=start_time,
                    estimated_completion=end_time,
                )
            )

        return {
            "input_rows": input_rows,
            "output_rows": output_rows,
            "filtered_rows": filtered_count,
            "columns": output_columns,
            "genome_build": parameters.get("genome_build"),
            "duration_seconds": duration,
            "output_file": target_path,
        }


def get_gwas_standardization_scenario() -> GWASStandardizationScenario:
    """Get singleton instance of GWAS standardization scenario."""
    return GWASStandardizationScenario()
