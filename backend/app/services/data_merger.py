"""
Data Merger Service

Provides functionality for merging multiple data files with different strategies:
- Vertical merge: Stack rows (append/concatenate)
- Horizontal merge: Join columns (by common key)
- Smart merge: Automatic detection and intelligent merging

Supports various file formats and handles edge cases like:
- Missing columns
- Duplicate keys
- Type mismatches
- Large files
"""

from typing import List, Dict, Optional, Any, Tuple
from pathlib import Path
from enum import Enum
import pandas as pd
import numpy as np
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


class MergeMode(str, Enum):
    """Data merge modes"""

    VERTICAL = "vertical"  # Stack rows (pd.concat axis=0)
    HORIZONTAL = "horizontal"  # Join columns (pd.merge)
    SMART = "smart"  # Auto-detect best strategy


class JoinType(str, Enum):
    """Join types for horizontal merge"""

    INNER = "inner"  # Only matching keys
    LEFT = "left"  # All from left, matching from right
    RIGHT = "right"  # All from right, matching from left
    OUTER = "outer"  # All keys from both


class DuplicateHandling(str, Enum):
    """How to handle duplicate keys/rows"""

    KEEP_FIRST = "keep_first"
    KEEP_LAST = "keep_last"
    KEEP_ALL = "keep_all"
    ERROR = "error"


@dataclass
class MergeConfig:
    """Configuration for data merge operation"""

    mode: MergeMode
    join_type: Optional[JoinType] = JoinType.OUTER
    key_columns: Optional[List[str]] = None  # For horizontal merge
    duplicate_handling: DuplicateHandling = DuplicateHandling.KEEP_ALL
    add_source_column: bool = True  # Add column indicating source file
    ignore_index: bool = False  # Reset index after merge
    sort_keys: bool = False  # Sort by key columns after merge
    fill_missing: Optional[Any] = None  # Value to fill missing data


@dataclass
class MergeResult:
    """Result of merge operation"""

    success: bool
    output_path: Optional[Path] = None
    n_rows: int = 0
    n_columns: int = 0
    n_files_merged: int = 0
    warnings: List[str] = None
    statistics: Dict[str, Any] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []
        if self.statistics is None:
            self.statistics = {}


class DataMerger:
    """
    Service for merging multiple data files.
    """

    def __init__(self):
        self.supported_formats = [".csv", ".tsv", ".txt", ".xlsx", ".parquet"]

    def merge_files(
        self,
        file_paths: List[Path],
        output_path: Path,
        config: MergeConfig,
    ) -> MergeResult:
        """
        Merge multiple files according to configuration.

        Args:
            file_paths: List of input file paths
            output_path: Path for merged output file
            config: Merge configuration

        Returns:
            MergeResult with success status and statistics
        """
        try:
            # Validate inputs
            validation_result = self._validate_inputs(file_paths, config)
            if not validation_result["valid"]:
                return MergeResult(success=False, warnings=validation_result["errors"])

            # Load all files
            dataframes = self._load_files(file_paths)

            # Add source tracking if requested
            if config.add_source_column:
                dataframes = self._add_source_tracking(dataframes, file_paths)

            # Perform merge based on mode
            if config.mode == MergeMode.VERTICAL:
                merged_df = self._vertical_merge(dataframes, config)
            elif config.mode == MergeMode.HORIZONTAL:
                merged_df = self._horizontal_merge(dataframes, config)
            elif config.mode == MergeMode.SMART:
                # Auto-detect best merge strategy
                config = self._detect_merge_strategy(dataframes, config)
                if config.mode == MergeMode.VERTICAL:
                    merged_df = self._vertical_merge(dataframes, config)
                else:
                    merged_df = self._horizontal_merge(dataframes, config)
            else:
                raise ValueError(f"Unknown merge mode: {config.mode}")

            # Post-processing
            merged_df = self._post_process(merged_df, config)

            # Save output
            self._save_output(merged_df, output_path)

            # Generate statistics
            stats = self._generate_statistics(dataframes, merged_df)

            return MergeResult(
                success=True,
                output_path=output_path,
                n_rows=len(merged_df),
                n_columns=len(merged_df.columns),
                n_files_merged=len(file_paths),
                statistics=stats,
            )

        except Exception as e:
            logger.error(f"Merge failed: {str(e)}")
            return MergeResult(success=False, warnings=[str(e)])

    def _validate_inputs(
        self, file_paths: List[Path], config: MergeConfig
    ) -> Dict[str, Any]:
        """Validate input files and configuration"""
        errors = []

        # Check minimum files
        if len(file_paths) < 2:
            errors.append("Need at least 2 files to merge")

        # Check files exist
        for path in file_paths:
            if not path.exists():
                errors.append(f"File not found: {path}")

        # Check file formats
        for path in file_paths:
            if path.suffix.lower() not in self.supported_formats:
                errors.append(f"Unsupported format: {path.suffix}")

        # Check horizontal merge has key columns
        if config.mode == MergeMode.HORIZONTAL and not config.key_columns:
            errors.append("Horizontal merge requires key_columns")

        return {"valid": len(errors) == 0, "errors": errors}

    def _load_files(self, file_paths: List[Path]) -> List[pd.DataFrame]:
        """Load all files into DataFrames"""
        dataframes = []

        for path in file_paths:
            try:
                df = self._load_single_file(path)
                dataframes.append(df)
                logger.info(f"Loaded {path}: {df.shape[0]} rows, {df.shape[1]} cols")
            except Exception as e:
                logger.error(f"Failed to load {path}: {e}")
                raise

        return dataframes

    def _load_single_file(self, path: Path) -> pd.DataFrame:
        """Load a single file into DataFrame"""
        suffix = path.suffix.lower()

        if suffix == ".csv":
            return pd.read_csv(path)
        elif suffix in [".tsv", ".txt"]:
            return pd.read_csv(path, sep="\t")
        elif suffix == ".xlsx":
            return pd.read_excel(path)
        elif suffix == ".parquet":
            return pd.read_parquet(path)
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

    def _add_source_tracking(
        self, dataframes: List[pd.DataFrame], file_paths: List[Path]
    ) -> List[pd.DataFrame]:
        """Add source file column to each DataFrame"""
        result = []

        for df, path in zip(dataframes, file_paths):
            df = df.copy()
            df["_source_file"] = path.name
            result.append(df)

        return result

    def _vertical_merge(
        self, dataframes: List[pd.DataFrame], config: MergeConfig
    ) -> pd.DataFrame:
        """
        Vertical merge: Stack rows (concatenate)

        Handles:
        - Missing columns (fill with NaN or config value)
        - Duplicate rows
        - Index reset
        """
        # Concatenate vertically
        merged = pd.concat(
            dataframes, axis=0, ignore_index=config.ignore_index, sort=False
        )

        # Handle duplicates
        if config.duplicate_handling == DuplicateHandling.KEEP_FIRST:
            merged = merged.drop_duplicates(keep="first")
        elif config.duplicate_handling == DuplicateHandling.KEEP_LAST:
            merged = merged.drop_duplicates(keep="last")
        elif config.duplicate_handling == DuplicateHandling.ERROR:
            if merged.duplicated().any():
                raise ValueError("Duplicate rows found")

        return merged

    def _horizontal_merge(
        self, dataframes: List[pd.DataFrame], config: MergeConfig
    ) -> pd.DataFrame:
        """
        Horizontal merge: Join by key columns

        Handles:
        - Multiple join types (inner/left/right/outer)
        - Multiple key columns
        - Missing keys
        - Duplicate keys
        """
        if not config.key_columns:
            raise ValueError("key_columns required for horizontal merge")

        # Start with first dataframe
        merged = dataframes[0]

        # Sequentially merge remaining dataframes
        for i, df in enumerate(dataframes[1:], start=1):
            try:
                merged = pd.merge(
                    merged,
                    df,
                    on=config.key_columns,
                    how=config.join_type.value,
                    suffixes=(f"_file{i}", f"_file{i+1}"),
                )
            except KeyError as e:
                logger.warning(f"Key column not found in file {i+1}: {e}")
                # Try to continue with available keys
                available_keys = [k for k in config.key_columns if k in df.columns]
                if available_keys:
                    merged = pd.merge(
                        merged,
                        df,
                        on=available_keys,
                        how=config.join_type.value,
                        suffixes=(f"_file{i}", f"_file{i+1}"),
                    )
                else:
                    raise ValueError(f"No common key columns found in file {i+1}")

        return merged

    def _detect_merge_strategy(
        self, dataframes: List[pd.DataFrame], config: MergeConfig
    ) -> MergeConfig:
        """
        Auto-detect best merge strategy (Smart merge)

        Heuristics:
        1. If all dataframes have same columns → Vertical merge
        2. If there are common columns that could be keys → Horizontal merge
        3. Otherwise → Vertical merge with column union
        """
        # Check if all have same columns
        first_cols = set(dataframes[0].columns)
        all_same = all(set(df.columns) == first_cols for df in dataframes[1:])

        if all_same:
            logger.info("Smart merge: Using VERTICAL (same columns)")
            config.mode = MergeMode.VERTICAL
            return config

        # Find common columns
        common_cols = first_cols.copy()
        for df in dataframes[1:]:
            common_cols &= set(df.columns)

        if len(common_cols) > 0:
            # Check if common columns could be good keys
            # (have unique or mostly unique values)
            potential_keys = []

            for col in common_cols:
                # Check uniqueness ratio in first dataframe
                n_unique = dataframes[0][col].nunique()
                n_total = len(dataframes[0])

                if n_unique / n_total > 0.8:  # >80% unique
                    potential_keys.append(col)

            if potential_keys:
                logger.info(f"Smart merge: Using HORIZONTAL with keys {potential_keys}")
                config.mode = MergeMode.HORIZONTAL
                config.key_columns = potential_keys[:1]  # Use most unique
                return config

        # Default to vertical merge
        logger.info("Smart merge: Using VERTICAL (column union)")
        config.mode = MergeMode.VERTICAL
        return config

    def _post_process(self, df: pd.DataFrame, config: MergeConfig) -> pd.DataFrame:
        """Post-processing steps"""

        # Fill missing values
        if config.fill_missing is not None:
            df = df.fillna(config.fill_missing)

        # Sort by key columns
        if config.sort_keys and config.key_columns:
            available_keys = [k for k in config.key_columns if k in df.columns]
            if available_keys:
                df = df.sort_values(available_keys)

        # Reset index if needed
        if config.ignore_index:
            df = df.reset_index(drop=True)

        return df

    def _save_output(self, df: pd.DataFrame, output_path: Path):
        """Save merged DataFrame to file"""
        suffix = output_path.suffix.lower()

        if suffix == ".csv":
            df.to_csv(output_path, index=False)
        elif suffix in [".tsv", ".txt"]:
            df.to_csv(output_path, sep="\t", index=False)
        elif suffix == ".xlsx":
            df.to_excel(output_path, index=False)
        elif suffix == ".parquet":
            df.to_parquet(output_path, index=False)
        else:
            # Default to CSV
            df.to_csv(output_path, index=False)

        logger.info(f"Saved merged data to {output_path}")

    def _generate_statistics(
        self, input_dfs: List[pd.DataFrame], merged_df: pd.DataFrame
    ) -> Dict[str, Any]:
        """Generate merge statistics"""

        input_stats = []
        for i, df in enumerate(input_dfs, 1):
            input_stats.append(
                {"file_index": i, "rows": len(df), "columns": len(df.columns)}
            )

        return {
            "input_files": input_stats,
            "output": {"rows": len(merged_df), "columns": len(merged_df.columns)},
            "total_input_rows": sum(len(df) for df in input_dfs),
            "rows_added": len(merged_df) - len(input_dfs[0]),
            "columns_added": len(merged_df.columns) - len(input_dfs[0].columns),
        }

    def preview_merge(
        self,
        file_paths: List[Path],
        config: MergeConfig,
        n_rows: int = 10,
    ) -> Dict[str, Any]:
        """
        Generate preview of merge operation without saving

        Args:
            file_paths: Input files
            config: Merge configuration
            n_rows: Number of rows to preview

        Returns:
            Dict with preview data and statistics
        """
        try:
            # Load sample from each file
            dataframes = []
            for path in file_paths:
                df = self._load_single_file(path)
                dataframes.append(df.head(n_rows * 2))  # Load extra for preview

            # Add source tracking if requested
            if config.add_source_column:
                dataframes = self._add_source_tracking(dataframes, file_paths)

            # Perform merge
            if config.mode == MergeMode.VERTICAL:
                merged_df = self._vertical_merge(dataframes, config)
            elif config.mode == MergeMode.HORIZONTAL:
                merged_df = self._horizontal_merge(dataframes, config)
            elif config.mode == MergeMode.SMART:
                config = self._detect_merge_strategy(dataframes, config)
                if config.mode == MergeMode.VERTICAL:
                    merged_df = self._vertical_merge(dataframes, config)
                else:
                    merged_df = self._horizontal_merge(dataframes, config)

            # Get preview
            preview_df = merged_df.head(n_rows)

            return {
                "success": True,
                "preview_data": {
                    "columns": list(preview_df.columns),
                    "rows": preview_df.values.tolist(),
                },
                "detected_mode": config.mode.value,
                "key_columns": config.key_columns,
                "estimated_rows": len(merged_df),
                "estimated_columns": len(merged_df.columns),
            }

        except Exception as e:
            logger.error(f"Preview failed: {e}")
            return {"success": False, "error": str(e)}


# Singleton instance
_data_merger = None


def get_data_merger() -> DataMerger:
    """Get singleton DataMerger instance"""
    global _data_merger
    if _data_merger is None:
        _data_merger = DataMerger()
    return _data_merger
