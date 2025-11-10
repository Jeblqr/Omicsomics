"""
Quick Visualizer Service

Automatically generates appropriate visualizations for data files
based on their content and structure. Supports smart detection of
data types and suggests best visualization methods.
"""

from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
from enum import Enum
import pandas as pd
import numpy as np
import logging
import json
import base64
from io import BytesIO

logger = logging.getLogger(__name__)

# Import visualization libraries
try:
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("Matplotlib/Seaborn not available")


class DataType(str, Enum):
    """Detected data types"""

    NUMERIC_MATRIX = "numeric_matrix"  # 数值矩阵
    TIME_SERIES = "time_series"  # 时间序列
    CATEGORICAL = "categorical"  # 分类数据
    MIXED = "mixed"  # 混合类型
    GENOMIC_INTERVALS = "genomic_intervals"  # 基因组区间
    EXPRESSION_MATRIX = "expression_matrix"  # 表达矩阵
    UNKNOWN = "unknown"


class ChartType(str, Enum):
    """Chart types"""

    HISTOGRAM = "histogram"
    SCATTER = "scatter"
    LINE = "line"
    BAR = "bar"
    BOXPLOT = "boxplot"
    HEATMAP = "heatmap"
    VIOLIN = "violin"
    CORRELATION = "correlation"


class QuickVisualizer:
    """Service for quick data visualization"""

    def __init__(self):
        self.supported_formats = [".csv", ".tsv", ".txt", ".xlsx"]
        self.max_rows_preview = 10000  # Max rows for quick viz

    def analyze_file(self, file_path: Path) -> Dict[str, Any]:
        """
        Analyze file and suggest visualizations

        Returns:
            Dict with data type, suggested charts, and column info
        """
        try:
            # Load data
            df = self._load_file(file_path)

            if df is None or df.empty:
                return {
                    "success": False,
                    "error": "Failed to load file or file is empty",
                }

            # Limit rows for performance
            if len(df) > self.max_rows_preview:
                df = df.head(self.max_rows_preview)
                limited = True
            else:
                limited = False

            # Detect data type
            data_type = self._detect_data_type(df)

            # Analyze columns
            column_info = self._analyze_columns(df)

            # Suggest visualizations
            suggestions = self._suggest_visualizations(df, data_type, column_info)

            return {
                "success": True,
                "data_type": data_type.value,
                "n_rows": len(df),
                "n_columns": len(df.columns),
                "limited": limited,
                "column_info": column_info,
                "suggestions": suggestions,
            }

        except Exception as e:
            logger.error(f"Failed to analyze file: {e}")
            return {"success": False, "error": str(e)}

    def generate_visualization(
        self, file_path: Path, chart_type: str, options: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Generate visualization

        Args:
            file_path: Path to data file
            chart_type: Type of chart to generate
            options: Visualization options (columns, limits, etc.)

        Returns:
            Dict with image data (base64) or error
        """
        if not MATPLOTLIB_AVAILABLE:
            return {"success": False, "error": "Visualization libraries not available"}

        try:
            # Load data
            df = self._load_file(file_path)

            if df is None or df.empty:
                return {"success": False, "error": "Failed to load file"}

            # Limit rows
            if len(df) > self.max_rows_preview:
                df = df.head(self.max_rows_preview)

            # Generate chart based on type
            chart_enum = ChartType(chart_type)
            options = options or {}

            if chart_enum == ChartType.HISTOGRAM:
                fig = self._create_histogram(df, options)
            elif chart_enum == ChartType.SCATTER:
                fig = self._create_scatter(df, options)
            elif chart_enum == ChartType.LINE:
                fig = self._create_line(df, options)
            elif chart_enum == ChartType.BAR:
                fig = self._create_bar(df, options)
            elif chart_enum == ChartType.BOXPLOT:
                fig = self._create_boxplot(df, options)
            elif chart_enum == ChartType.HEATMAP:
                fig = self._create_heatmap(df, options)
            elif chart_enum == ChartType.VIOLIN:
                fig = self._create_violin(df, options)
            elif chart_enum == ChartType.CORRELATION:
                fig = self._create_correlation(df, options)
            else:
                return {
                    "success": False,
                    "error": f"Unsupported chart type: {chart_type}",
                }

            # Convert to base64
            buffer = BytesIO()
            fig.savefig(buffer, format="png", dpi=100, bbox_inches="tight")
            buffer.seek(0)
            image_base64 = base64.b64encode(buffer.read()).decode()
            plt.close(fig)

            return {
                "success": True,
                "image": f"data:image/png;base64,{image_base64}",
                "chart_type": chart_type,
            }

        except Exception as e:
            logger.error(f"Failed to generate visualization: {e}")
            return {"success": False, "error": str(e)}

    def _load_file(self, file_path: Path) -> Optional[pd.DataFrame]:
        """Load file into DataFrame"""
        suffix = file_path.suffix.lower()

        try:
            if suffix == ".csv":
                return pd.read_csv(file_path)
            elif suffix in [".tsv", ".txt"]:
                return pd.read_csv(file_path, sep="\t")
            elif suffix == ".xlsx":
                return pd.read_excel(file_path)
            else:
                return None
        except Exception as e:
            logger.error(f"Failed to load {file_path}: {e}")
            return None

    def _detect_data_type(self, df: pd.DataFrame) -> DataType:
        """Detect overall data type"""

        # Check for expression matrix pattern (many numeric columns)
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 10 and len(numeric_cols) / len(df.columns) > 0.8:
            return DataType.EXPRESSION_MATRIX

        # Check for time series
        date_cols = df.select_dtypes(include=["datetime64"]).columns
        if len(date_cols) > 0:
            return DataType.TIME_SERIES

        # Check column names for time patterns
        time_keywords = ["time", "date", "year", "month", "day"]
        if any(
            any(kw in str(col).lower() for kw in time_keywords) for col in df.columns
        ):
            return DataType.TIME_SERIES

        # Check for genomic intervals
        genomic_keywords = ["chr", "chromosome", "start", "end", "pos", "position"]
        if (
            sum(
                any(kw in str(col).lower() for kw in genomic_keywords)
                for col in df.columns
            )
            >= 2
        ):
            return DataType.GENOMIC_INTERVALS

        # Check for mostly numeric
        if len(numeric_cols) / len(df.columns) > 0.7:
            return DataType.NUMERIC_MATRIX

        # Check for categorical
        categorical_cols = df.select_dtypes(include=["object", "category"]).columns
        if len(categorical_cols) / len(df.columns) > 0.5:
            return DataType.CATEGORICAL

        # Mixed
        if len(numeric_cols) > 0 and len(categorical_cols) > 0:
            return DataType.MIXED

        return DataType.UNKNOWN

    def _analyze_columns(self, df: pd.DataFrame) -> List[Dict[str, Any]]:
        """Analyze each column"""
        column_info = []

        for col in df.columns:
            col_data = df[col]

            info = {
                "name": str(col),
                "dtype": str(col_data.dtype),
                "n_unique": int(col_data.nunique()),
                "n_missing": int(col_data.isnull().sum()),
                "is_numeric": pd.api.types.is_numeric_dtype(col_data),
                "is_categorical": pd.api.types.is_categorical_dtype(col_data)
                or col_data.dtype == "object",
            }

            if info["is_numeric"]:
                info["min"] = float(col_data.min()) if not col_data.empty else None
                info["max"] = float(col_data.max()) if not col_data.empty else None
                info["mean"] = float(col_data.mean()) if not col_data.empty else None
                info["median"] = (
                    float(col_data.median()) if not col_data.empty else None
                )

            column_info.append(info)

        return column_info

    def _suggest_visualizations(
        self, df: pd.DataFrame, data_type: DataType, column_info: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Suggest appropriate visualizations"""
        suggestions = []

        numeric_cols = [c["name"] for c in column_info if c["is_numeric"]]
        categorical_cols = [c["name"] for c in column_info if c["is_categorical"]]

        # Expression matrix / Numeric matrix
        if data_type in [DataType.EXPRESSION_MATRIX, DataType.NUMERIC_MATRIX]:
            suggestions.append(
                {
                    "chart_type": ChartType.HEATMAP.value,
                    "title": "Heatmap",
                    "description": "Visualize numeric data as heatmap",
                    "options": {},
                    "priority": 1,
                }
            )

            suggestions.append(
                {
                    "chart_type": ChartType.CORRELATION.value,
                    "title": "Correlation Matrix",
                    "description": "Show correlations between columns",
                    "options": {},
                    "priority": 2,
                }
            )

            if len(numeric_cols) > 0:
                suggestions.append(
                    {
                        "chart_type": ChartType.BOXPLOT.value,
                        "title": "Box Plot",
                        "description": "Distribution of numeric columns",
                        "options": {"columns": numeric_cols[:10]},
                        "priority": 3,
                    }
                )

        # Time series
        elif data_type == DataType.TIME_SERIES:
            suggestions.append(
                {
                    "chart_type": ChartType.LINE.value,
                    "title": "Line Chart",
                    "description": "Visualize trends over time",
                    "options": {},
                    "priority": 1,
                }
            )

        # Categorical
        elif data_type == DataType.CATEGORICAL:
            if len(categorical_cols) > 0:
                suggestions.append(
                    {
                        "chart_type": ChartType.BAR.value,
                        "title": "Bar Chart",
                        "description": "Count by category",
                        "options": {"column": categorical_cols[0]},
                        "priority": 1,
                    }
                )

        # Mixed or general suggestions
        if len(numeric_cols) >= 2:
            suggestions.append(
                {
                    "chart_type": ChartType.SCATTER.value,
                    "title": "Scatter Plot",
                    "description": f"Plot {numeric_cols[0]} vs {numeric_cols[1]}",
                    "options": {
                        "x_column": numeric_cols[0],
                        "y_column": numeric_cols[1],
                    },
                    "priority": 2,
                }
            )

        if len(numeric_cols) >= 1:
            suggestions.append(
                {
                    "chart_type": ChartType.HISTOGRAM.value,
                    "title": "Histogram",
                    "description": f"Distribution of {numeric_cols[0]}",
                    "options": {"column": numeric_cols[0]},
                    "priority": 3,
                }
            )

        # Sort by priority
        suggestions.sort(key=lambda x: x["priority"])

        return suggestions

    # Chart creation methods
    def _create_histogram(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create histogram"""
        fig, ax = plt.subplots(figsize=(10, 6))

        column = options.get("column")
        if not column:
            # Use first numeric column
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) == 0:
                raise ValueError("No numeric columns found")
            column = numeric_cols[0]

        data = df[column].dropna()
        ax.hist(data, bins=options.get("bins", 30), edgecolor="black", alpha=0.7)
        ax.set_xlabel(column)
        ax.set_ylabel("Frequency")
        ax.set_title(f"Distribution of {column}")
        ax.grid(True, alpha=0.3)

        return fig

    def _create_scatter(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create scatter plot"""
        fig, ax = plt.subplots(figsize=(10, 6))

        x_col = options.get("x_column")
        y_col = options.get("y_column")

        if not x_col or not y_col:
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) < 2:
                raise ValueError("Need at least 2 numeric columns")
            x_col = numeric_cols[0]
            y_col = numeric_cols[1]

        ax.scatter(df[x_col], df[y_col], alpha=0.6)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_title(f"{y_col} vs {x_col}")
        ax.grid(True, alpha=0.3)

        return fig

    def _create_line(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create line chart"""
        fig, ax = plt.subplots(figsize=(12, 6))

        numeric_cols = df.select_dtypes(include=[np.number]).columns
        y_cols = options.get("columns", list(numeric_cols[:5]))  # Max 5 lines

        for col in y_cols:
            if col in df.columns:
                ax.plot(df.index, df[col], label=col, marker="o", markersize=3)

        ax.set_xlabel("Index")
        ax.set_ylabel("Value")
        ax.set_title("Line Chart")
        ax.legend()
        ax.grid(True, alpha=0.3)

        return fig

    def _create_bar(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create bar chart"""
        fig, ax = plt.subplots(figsize=(12, 6))

        column = options.get("column")
        if not column:
            # Use first categorical column
            cat_cols = df.select_dtypes(include=["object", "category"]).columns
            if len(cat_cols) == 0:
                raise ValueError("No categorical columns found")
            column = cat_cols[0]

        value_counts = df[column].value_counts().head(20)  # Top 20
        ax.bar(range(len(value_counts)), value_counts.values)
        ax.set_xticks(range(len(value_counts)))
        ax.set_xticklabels(value_counts.index, rotation=45, ha="right")
        ax.set_ylabel("Count")
        ax.set_title(f"Count by {column}")
        ax.grid(True, axis="y", alpha=0.3)

        return fig

    def _create_boxplot(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create box plot"""
        fig, ax = plt.subplots(figsize=(12, 6))

        numeric_cols = df.select_dtypes(include=[np.number]).columns
        columns = options.get("columns", list(numeric_cols[:10]))  # Max 10

        data_to_plot = [df[col].dropna() for col in columns if col in df.columns]
        ax.boxplot(data_to_plot, labels=columns)
        ax.set_ylabel("Value")
        ax.set_title("Box Plot")
        ax.grid(True, axis="y", alpha=0.3)
        plt.xticks(rotation=45, ha="right")

        return fig

    def _create_heatmap(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create heatmap"""
        fig, ax = plt.subplots(figsize=(12, 10))

        # Use only numeric columns
        numeric_df = df.select_dtypes(include=[np.number])

        # Limit size for performance
        if len(numeric_df) > 100:
            numeric_df = numeric_df.head(100)
        if len(numeric_df.columns) > 50:
            numeric_df = numeric_df.iloc[:, :50]

        sns.heatmap(numeric_df, cmap="viridis", ax=ax, cbar_kws={"label": "Value"})
        ax.set_title("Heatmap")

        return fig

    def _create_violin(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create violin plot"""
        fig, ax = plt.subplots(figsize=(12, 6))

        numeric_cols = df.select_dtypes(include=[np.number]).columns
        columns = options.get("columns", list(numeric_cols[:5]))

        data_to_plot = [df[col].dropna() for col in columns if col in df.columns]
        ax.violinplot(data_to_plot, positions=range(len(columns)), showmeans=True)
        ax.set_xticks(range(len(columns)))
        ax.set_xticklabels(columns, rotation=45, ha="right")
        ax.set_ylabel("Value")
        ax.set_title("Violin Plot")
        ax.grid(True, axis="y", alpha=0.3)

        return fig

    def _create_correlation(self, df: pd.DataFrame, options: Dict) -> plt.Figure:
        """Create correlation matrix"""
        fig, ax = plt.subplots(figsize=(12, 10))

        numeric_df = df.select_dtypes(include=[np.number])

        # Limit columns
        if len(numeric_df.columns) > 20:
            numeric_df = numeric_df.iloc[:, :20]

        corr_matrix = numeric_df.corr()

        sns.heatmap(
            corr_matrix,
            annot=True,
            fmt=".2f",
            cmap="coolwarm",
            center=0,
            ax=ax,
            square=True,
            cbar_kws={"label": "Correlation"},
        )
        ax.set_title("Correlation Matrix")

        return fig


# Singleton instance
_quick_visualizer = None


def get_quick_visualizer() -> QuickVisualizer:
    """Get singleton QuickVisualizer instance"""
    global _quick_visualizer
    if _quick_visualizer is None:
        _quick_visualizer = QuickVisualizer()
    return _quick_visualizer
