"""
Interactive format converter base framework.

Provides base classes and utilities for user-guided format conversions
that require parameter configuration, column mapping, or data interpretation.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from enum import Enum
from dataclasses import dataclass, field
from datetime import datetime
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class ParameterType(Enum):
    """Types of conversion parameters."""

    TEXT = "text"
    NUMBER = "number"
    SELECT = "select"
    MULTI_SELECT = "multi_select"
    BOOLEAN = "boolean"
    FILE = "file"
    COLUMN_MAPPING = "column_mapping"
    THRESHOLD = "threshold"


class ValidationLevel(Enum):
    """Validation message severity levels."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass
class ConversionParameter:
    """
    Definition of a conversion parameter.

    Attributes:
        name: Parameter identifier
        label: Human-readable label
        type: Parameter type (from ParameterType enum)
        required: Whether parameter is required
        default: Default value
        options: List of options for select/multi_select types
        description: Help text for the parameter
        validation: Validation rules
    """

    name: str
    label: str
    type: ParameterType
    required: bool = False
    default: Any = None
    options: Optional[List[Dict[str, Any]]] = None
    description: str = ""
    validation: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API responses."""
        return {
            "name": self.name,
            "label": self.label,
            "type": self.type.value,
            "required": self.required,
            "default": self.default,
            "options": self.options,
            "description": self.description,
            "validation": self.validation,
        }


@dataclass
class ValidationMessage:
    """
    A validation message for user feedback.

    Attributes:
        level: Severity level
        message: Human-readable message
        field: Related parameter field (optional)
        details: Additional details
    """

    level: ValidationLevel
    message: str
    field: Optional[str] = None
    details: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API responses."""
        return {
            "level": self.level.value,
            "message": self.message,
            "field": self.field,
            "details": self.details,
        }


@dataclass
class ConversionPreview:
    """
    Preview of conversion results.

    Attributes:
        input_info: Information about input file
        output_info: Information about expected output
        sample_data: Sample rows of output data
        validation_messages: List of validation messages
        estimated_time: Estimated conversion time in seconds
    """

    input_info: Dict[str, Any]
    output_info: Dict[str, Any]
    sample_data: Optional[List[Dict[str, Any]]] = None
    validation_messages: List[ValidationMessage] = field(default_factory=list)
    estimated_time: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API responses."""
        return {
            "input_info": self.input_info,
            "output_info": self.output_info,
            "sample_data": self.sample_data,
            "validation_messages": [msg.to_dict() for msg in self.validation_messages],
            "estimated_time": self.estimated_time,
        }


@dataclass
class ConversionProgress:
    """
    Progress tracking for conversion.

    Attributes:
        stage: Current conversion stage
        progress: Progress percentage (0-100)
        message: Current status message
        started_at: Start timestamp
        estimated_completion: Estimated completion timestamp
    """

    stage: str
    progress: float
    message: str
    started_at: datetime
    estimated_completion: Optional[datetime] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API responses."""
        return {
            "stage": self.stage,
            "progress": self.progress,
            "message": self.message,
            "started_at": self.started_at.isoformat(),
            "estimated_completion": (
                self.estimated_completion.isoformat()
                if self.estimated_completion
                else None
            ),
        }


class ConversionScenario(ABC):
    """
    Base class for interactive conversion scenarios.

    Each scenario represents a specific conversion workflow that requires
    user interaction (e.g., GWAS standardization, single-cell import).
    """

    def __init__(self):
        self.name = self.__class__.__name__
        self.parameters: List[ConversionParameter] = []
        self._initialize_parameters()

    @abstractmethod
    def _initialize_parameters(self):
        """Initialize the parameters required for this conversion scenario."""
        pass

    @abstractmethod
    def detect_format(self, file_path: str) -> bool:
        """
        Detect if the input file matches this conversion scenario.

        Args:
            file_path: Path to input file

        Returns:
            True if file matches this scenario
        """
        pass

    @abstractmethod
    def validate_input(self, file_path: str) -> List[ValidationMessage]:
        """
        Validate the input file.

        Args:
            file_path: Path to input file

        Returns:
            List of validation messages
        """
        pass

    @abstractmethod
    def validate_parameters(
        self, parameters: Dict[str, Any]
    ) -> List[ValidationMessage]:
        """
        Validate user-provided parameters.

        Args:
            parameters: Dictionary of parameter values

        Returns:
            List of validation messages
        """
        pass

    @abstractmethod
    def generate_preview(
        self, file_path: str, parameters: Dict[str, Any]
    ) -> ConversionPreview:
        """
        Generate a preview of the conversion results.

        Args:
            file_path: Path to input file
            parameters: Dictionary of parameter values

        Returns:
            ConversionPreview object
        """
        pass

    @abstractmethod
    def convert(
        self,
        source_path: str,
        target_path: str,
        parameters: Dict[str, Any],
        progress_callback: Optional[callable] = None,
    ) -> Dict[str, Any]:
        """
        Perform the actual conversion.

        Args:
            source_path: Path to input file
            target_path: Path to output file
            parameters: Dictionary of parameter values
            progress_callback: Optional callback for progress updates

        Returns:
            Dictionary with conversion results and statistics
        """
        pass

    def get_parameter_definitions(self) -> List[Dict[str, Any]]:
        """
        Get parameter definitions for API responses.

        Returns:
            List of parameter definitions as dictionaries
        """
        return [param.to_dict() for param in self.parameters]

    def get_scenario_info(self) -> Dict[str, Any]:
        """
        Get scenario information for API responses.

        Returns:
            Dictionary with scenario metadata
        """
        return {
            "name": self.name,
            "description": self.__doc__ or "",
            "parameters": self.get_parameter_definitions(),
        }


class InteractiveConverter:
    """
    Main interactive converter class.

    Manages multiple conversion scenarios and routes conversions
    to the appropriate scenario handler.
    """

    def __init__(self):
        self.scenarios: Dict[str, ConversionScenario] = {}
        self._initialize_scenarios()

    def _initialize_scenarios(self):
        """Initialize built-in scenarios"""
        from .scenarios import (
            get_gwas_standardization_scenario,
            get_expression_matrix_standardization_scenario,
            get_single_cell_import_scenario,
            get_proteomics_standardization_scenario,
            get_metabolomics_standardization_scenario,
            get_epigenomics_conversion_scenario,
            get_multiomics_integration_scenario,
            get_network_pathway_formatting_scenario,
            get_clinical_data_standardization_scenario,
        )

        # Register all interactive scenarios
        self.register_scenario(get_gwas_standardization_scenario())
        self.register_scenario(get_expression_matrix_standardization_scenario())
        self.register_scenario(get_single_cell_import_scenario())
        self.register_scenario(get_proteomics_standardization_scenario())
        self.register_scenario(get_metabolomics_standardization_scenario())
        self.register_scenario(get_epigenomics_conversion_scenario())
        self.register_scenario(get_multiomics_integration_scenario())
        self.register_scenario(get_network_pathway_formatting_scenario())
        self.register_scenario(get_clinical_data_standardization_scenario())

    def register_scenario(self, scenario_id: str, scenario: ConversionScenario):
        """
        Register a new conversion scenario.

        Args:
            scenario_id: Unique identifier for the scenario
            scenario: ConversionScenario instance
        """
        self.scenarios[scenario_id] = scenario
        logger.info(f"Registered conversion scenario: {scenario_id}")

    def detect_scenario(self, file_path: str) -> Optional[str]:
        """
        Detect which conversion scenario matches the input file.

        Args:
            file_path: Path to input file

        Returns:
            Scenario ID if match found, None otherwise
        """
        for scenario_id, scenario in self.scenarios.items():
            try:
                if scenario.detect_format(file_path):
                    return scenario_id
            except Exception as e:
                logger.warning(f"Error detecting scenario {scenario_id}: {e}")
                continue

        return None

    def get_scenario(self, scenario_id: str) -> Optional[ConversionScenario]:
        """
        Get a conversion scenario by ID.

        Args:
            scenario_id: Scenario identifier

        Returns:
            ConversionScenario instance or None
        """
        return self.scenarios.get(scenario_id)

    def list_scenarios(self) -> List[Dict[str, Any]]:
        """
        List all available conversion scenarios.

        Returns:
            List of scenario information dictionaries
        """
        return [
            {"id": scenario_id, **scenario.get_scenario_info()}
            for scenario_id, scenario in self.scenarios.items()
        ]

    def validate_file(
        self, file_path: str, scenario_id: str
    ) -> List[ValidationMessage]:
        """
        Validate an input file for a specific scenario.

        Args:
            file_path: Path to input file
            scenario_id: Scenario identifier

        Returns:
            List of validation messages
        """
        scenario = self.get_scenario(scenario_id)
        if not scenario:
            return [
                ValidationMessage(
                    level=ValidationLevel.ERROR,
                    message=f"Unknown scenario: {scenario_id}",
                )
            ]

        return scenario.validate_input(file_path)

    def preview_conversion(
        self, file_path: str, scenario_id: str, parameters: Dict[str, Any]
    ) -> ConversionPreview:
        """
        Generate a preview of conversion results.

        Args:
            file_path: Path to input file
            scenario_id: Scenario identifier
            parameters: Conversion parameters

        Returns:
            ConversionPreview object
        """
        scenario = self.get_scenario(scenario_id)
        if not scenario:
            raise ValueError(f"Unknown scenario: {scenario_id}")

        # Validate parameters first
        validation_messages = scenario.validate_parameters(parameters)

        # Generate preview
        preview = scenario.generate_preview(file_path, parameters)

        # Add parameter validation messages
        preview.validation_messages.extend(validation_messages)

        return preview

    def convert(
        self,
        source_path: str,
        target_path: str,
        scenario_id: str,
        parameters: Dict[str, Any],
        progress_callback: Optional[callable] = None,
    ) -> Dict[str, Any]:
        """
        Perform interactive conversion.

        Args:
            source_path: Path to input file
            target_path: Path to output file
            scenario_id: Scenario identifier
            parameters: Conversion parameters
            progress_callback: Optional callback for progress updates

        Returns:
            Dictionary with conversion results
        """
        scenario = self.get_scenario(scenario_id)
        if not scenario:
            raise ValueError(f"Unknown scenario: {scenario_id}")

        # Validate parameters
        validation_messages = scenario.validate_parameters(parameters)
        errors = [
            msg for msg in validation_messages if msg.level == ValidationLevel.ERROR
        ]

        if errors:
            raise ValueError(f"Parameter validation failed: {errors[0].message}")

        # Perform conversion
        result = scenario.convert(
            source_path, target_path, parameters, progress_callback
        )

        return result


# Singleton instance
_interactive_converter_instance = None


def get_interactive_converter() -> InteractiveConverter:
    """Get singleton instance of InteractiveConverter."""
    global _interactive_converter_instance
    if _interactive_converter_instance is None:
        _interactive_converter_instance = InteractiveConverter()
    return _interactive_converter_instance


# Utility functions for common operations


def detect_delimiter(file_path: str, max_lines: int = 10) -> str:
    """
    Detect the delimiter used in a text file.

    Args:
        file_path: Path to file
        max_lines: Maximum lines to check

    Returns:
        Detected delimiter character
    """
    common_delimiters = [",", "\t", " ", "|", ";"]
    delimiter_counts = {d: 0 for d in common_delimiters}

    with open(file_path, "r") as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            for delimiter in common_delimiters:
                delimiter_counts[delimiter] += line.count(delimiter)

    # Return delimiter with highest count
    return max(delimiter_counts.items(), key=lambda x: x[1])[0]


def preview_dataframe(df: pd.DataFrame, max_rows: int = 10) -> List[Dict[str, Any]]:
    """
    Convert DataFrame to list of dictionaries for preview.

    Args:
        df: Input DataFrame
        max_rows: Maximum rows to include

    Returns:
        List of row dictionaries
    """
    preview_df = df.head(max_rows)
    return preview_df.to_dict("records")


def infer_column_types(df: pd.DataFrame) -> Dict[str, str]:
    """
    Infer semantic types of DataFrame columns.

    Args:
        df: Input DataFrame

    Returns:
        Dictionary mapping column names to semantic types
    """
    column_types = {}

    for col in df.columns:
        dtype = str(df[col].dtype)

        # Check for numeric types
        if "int" in dtype or "float" in dtype:
            column_types[col] = "numeric"
        # Check for datetime
        elif "datetime" in dtype:
            column_types[col] = "datetime"
        # Check for boolean
        elif dtype == "bool":
            column_types[col] = "boolean"
        # Default to text
        else:
            column_types[col] = "text"

    return column_types
