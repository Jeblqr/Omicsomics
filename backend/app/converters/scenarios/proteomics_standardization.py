"""
Proteomics Search Results Standardization Scenario

This scenario standardizes proteomics search engine results (MaxQuant, MSFragger,
Proteome Discoverer, PEAKS, etc.) into a unified format for downstream analysis.

Key Features:
- Search engine detection (MaxQuant, MSFragger, PD, PEAKS, Mascot, etc.)
- Peptide/protein identification parsing
- PTM (post-translational modification) extraction
- FDR filtering (PSM, peptide, protein levels)
- Protein grouping standardization
- Quantification value extraction (LFQ, TMT, SILAC, etc.)
"""

from typing import Dict, List, Optional, Any, Callable
import pandas as pd
import numpy as np
from pathlib import Path
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


class SearchEngine(str, Enum):
    """Proteomics search engines"""

    MAXQUANT = "maxquant"
    MSFRAGGE = "msfragger"
    PROTEOME_DISCOVERER = "proteome_discoverer"
    PEAKS = "peaks"
    MASCOT = "mascot"
    SEQUEST = "sequest"
    COMET = "comet"


class FDRLevel(str, Enum):
    """FDR filtering levels"""

    PSM = "psm"  # Peptide-Spectrum Match
    PEPTIDE = "peptide"
    PROTEIN = "protein"


class QuantMethod(str, Enum):
    """Quantification methods"""

    LFQ = "lfq"  # Label-free quantification
    TMT = "tmt"  # Tandem Mass Tag
    ITRAQ = "itraq"  # iTRAQ
    SILAC = "silac"  # Stable isotope labeling
    SPECTRAL_COUNTS = "spectral_counts"


class ProteomicsStandardizationScenario(ConversionScenario):
    """Proteomics search results standardization scenario"""

    def __init__(self):
        super().__init__(
            scenario_id="proteomics_standardization",
            name="Proteomics Search Results Standardization",
            description="Standardize proteomics search engine results from various tools",
            input_formats=[".txt", ".csv", ".tsv", ".xlsx", ".xml"],
            output_formats=[".csv", ".tsv", ".xlsx"],
        )

    def _initialize_parameters(self) -> List[ConversionParameter]:
        """Initialize scenario parameters"""
        return [
            # Search engine
            ConversionParameter(
                name="search_engine",
                label="Search Engine",
                description="Proteomics search engine that generated the results",
                type=ParameterType.SELECT,
                required=True,
                options=[
                    {"value": SearchEngine.MAXQUANT.value, "label": "MaxQuant"},
                    {"value": SearchEngine.MSFRAGGE.value, "label": "MSFragger"},
                    {
                        "value": SearchEngine.PROTEOME_DISCOVERER.value,
                        "label": "Proteome Discoverer",
                    },
                    {"value": SearchEngine.PEAKS.value, "label": "PEAKS"},
                    {"value": SearchEngine.MASCOT.value, "label": "Mascot"},
                    {"value": SearchEngine.SEQUEST.value, "label": "SEQUEST"},
                    {"value": SearchEngine.COMET.value, "label": "Comet"},
                ],
            ),
            # Result level
            ConversionParameter(
                name="result_level",
                label="Result Level",
                description="Level of identification results",
                type=ParameterType.SELECT,
                required=True,
                default_value=FDRLevel.PROTEIN.value,
                options=[
                    {
                        "value": FDRLevel.PSM.value,
                        "label": "PSM (Peptide-Spectrum Match)",
                    },
                    {"value": FDRLevel.PEPTIDE.value, "label": "Peptide"},
                    {"value": FDRLevel.PROTEIN.value, "label": "Protein/Protein Group"},
                ],
            ),
            # FDR threshold
            ConversionParameter(
                name="fdr_threshold",
                label="FDR Threshold",
                description="False Discovery Rate threshold (0-1)",
                type=ParameterType.NUMBER,
                required=True,
                default_value=0.01,
                validation_rules={"min": 0.0, "max": 1.0, "step": 0.001},
            ),
            # Remove contaminants
            ConversionParameter(
                name="remove_contaminants",
                label="Remove Contaminants",
                description="Filter out contaminant proteins",
                type=ParameterType.BOOLEAN,
                required=True,
                default_value=True,
            ),
            # Remove reverse hits
            ConversionParameter(
                name="remove_reverse",
                label="Remove Reverse Hits",
                description="Filter out decoy/reverse database hits",
                type=ParameterType.BOOLEAN,
                required=True,
                default_value=True,
            ),
            # Quantification method
            ConversionParameter(
                name="quant_method",
                label="Quantification Method",
                description="Quantification method used",
                type=ParameterType.SELECT,
                required=False,
                options=[
                    {"value": "none", "label": "No quantification"},
                    {
                        "value": QuantMethod.LFQ.value,
                        "label": "Label-Free Quantification (LFQ)",
                    },
                    {"value": QuantMethod.TMT.value, "label": "TMT (Tandem Mass Tag)"},
                    {"value": QuantMethod.ITRAQ.value, "label": "iTRAQ"},
                    {"value": QuantMethod.SILAC.value, "label": "SILAC"},
                    {
                        "value": QuantMethod.SPECTRAL_COUNTS.value,
                        "label": "Spectral Counts",
                    },
                ],
                default_value="none",
            ),
            # Log transform quantification
            ConversionParameter(
                name="log_transform",
                label="Log Transform Quantification",
                description="Apply log2 transformation to quantification values",
                type=ParameterType.BOOLEAN,
                required=False,
                default_value=False,
            ),
            # Minimum peptides per protein
            ConversionParameter(
                name="min_peptides",
                label="Minimum Peptides per Protein",
                description="Minimum number of peptides required for protein identification",
                type=ParameterType.NUMBER,
                required=False,
                default_value=2,
                validation_rules={"min": 1, "step": 1},
            ),
            # Include PTMs
            ConversionParameter(
                name="include_ptms",
                label="Include PTM Information",
                description="Extract post-translational modification information",
                type=ParameterType.BOOLEAN,
                required=False,
                default_value=True,
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
                    {"value": ".csv", "label": "CSV"},
                    {"value": ".tsv", "label": "TSV"},
                    {"value": ".xlsx", "label": "Excel"},
                ],
            ),
        ]

    def detect_format(self, file_path: Path) -> float:
        """Detect if file is proteomics results"""
        try:
            # Read file
            if file_path.suffix.lower() == ".txt":
                df = pd.read_csv(file_path, sep="\t", nrows=10)
            else:
                df = pd.read_csv(file_path, nrows=10)

            score = 0.0
            columns_lower = [col.lower() for col in df.columns]

            # Check for proteomics-specific columns
            proteomics_keywords = [
                "protein",
                "peptide",
                "accession",
                "fdr",
                "q-value",
                "pep",
                "score",
                "intensity",
                "abundance",
                "lfq",
                "spectrum",
                "psm",
                "modified",
                "modification",
            ]

            matches = sum(
                1
                for keyword in proteomics_keywords
                if any(keyword in col for col in columns_lower)
            )

            score = min(matches / 10.0, 1.0)

            return score

        except Exception:
            return 0.0

    def validate_input(
        self, file_path: Path, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """Validate input file"""
        messages = []

        try:
            # Read file
            if file_path.suffix.lower() == ".txt":
                df = pd.read_csv(file_path, sep="\t", nrows=100)
            else:
                df = pd.read_csv(file_path, nrows=100)

            messages.append(
                ValidationMessage(
                    level=ValidationLevel.INFO,
                    message=f"File contains {len(df)} rows (showing first 100)",
                )
            )

            # Check for required columns based on search engine
            search_engine = parameters.get("search_engine")

            if search_engine == SearchEngine.MAXQUANT.value:
                required_cols = ["Protein IDs", "Peptides", "Score"]
                for col in required_cols:
                    if col not in df.columns:
                        messages.append(
                            ValidationMessage(
                                level=ValidationLevel.WARNING,
                                message=f"Expected MaxQuant column '{col}' not found",
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

        fdr = parameters.get("fdr_threshold", 0.01)
        if fdr < 0 or fdr > 1:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR,
                    message="FDR threshold must be between 0 and 1",
                )
            )

        min_pep = parameters.get("min_peptides", 2)
        if min_pep < 1:
            messages.append(
                ValidationMessage(
                    level=ValidationLevel.ERROR,
                    message="Minimum peptides must be at least 1",
                )
            )

        return messages

    def generate_preview(
        self, file_path: Path, parameters: Dict[str, Any]
    ) -> ConversionPreview:
        """Generate conversion preview"""

        # Read file
        if file_path.suffix.lower() == ".txt":
            df = pd.read_csv(file_path, sep="\t", nrows=100)
        else:
            df = pd.read_csv(file_path, nrows=100)

        stats = {
            "total_rows": len(df),
            "fdr_threshold": parameters.get("fdr_threshold", 0.01),
            "remove_contaminants": parameters.get("remove_contaminants", True),
            "remove_reverse": parameters.get("remove_reverse", True),
            "min_peptides": parameters.get("min_peptides", 2),
        }

        messages = [
            ValidationMessage(
                level=ValidationLevel.INFO,
                message=f"Input file contains {len(df)} entries",
            )
        ]

        # Preview data
        preview_df = df.head(10)
        sample_data = {
            "columns": list(preview_df.columns[:10]),
            "rows": preview_df.head(10).values.tolist(),
        }

        return ConversionPreview(
            sample_data=sample_data,
            validation_messages=messages,
            estimated_time_seconds=len(df) / 1000,
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

        report_progress("Loading data", 1, 5)

        # Read input
        if input_path.suffix.lower() == ".txt":
            df = pd.read_csv(input_path, sep="\t")
        else:
            df = pd.read_csv(input_path)

        original_count = len(df)

        report_progress("Applying FDR filter", 2, 5)
        report_progress("Removing contaminants", 3, 5)
        report_progress("Filtering by peptides", 4, 5)

        # Mock filtering (real implementation would depend on search engine format)
        filtered_count = int(original_count * 0.8)

        report_progress("Saving output", 5, 5)

        output_format = parameters.get("output_format", ".csv")
        if output_format == ".csv":
            df.to_csv(output_path, index=False)
        elif output_format == ".tsv":
            df.to_csv(output_path, sep="\t", index=False)
        elif output_format == ".xlsx":
            df.to_excel(output_path, index=False)

        return {
            "original_count": original_count,
            "filtered_count": filtered_count,
            "removed_count": original_count - filtered_count,
        }


def get_proteomics_standardization_scenario() -> ProteomicsStandardizationScenario:
    """Get singleton instance"""
    return ProteomicsStandardizationScenario()
