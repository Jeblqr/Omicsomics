#!/usr/bin/env python3
"""
Quick Visualization Test Script

Tests the quick visualization API endpoints with various data types.
"""

import requests
import json
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import base64
from io import BytesIO

API_BASE_URL = "http://localhost:8000"


def create_test_data():
    """Create various test datasets"""
    test_dir = Path(tempfile.mkdtemp(prefix="quick_viz_test_"))
    print(f"📁 Creating test data in: {test_dir}")

    # 1. Expression matrix
    np.random.seed(42)
    genes = [f"Gene_{i}" for i in range(100)]
    samples = [f"Sample_{i}" for i in range(20)]
    expr_data = pd.DataFrame(
        np.random.lognormal(2, 1, (100, 20)), index=genes, columns=samples
    )
    expr_file = test_dir / "expression_matrix.csv"
    expr_data.to_csv(expr_file)
    print(f"  ✓ Created expression matrix: {expr_file}")

    # 2. Time series
    dates = pd.date_range("2024-01-01", periods=100, freq="D")
    ts_data = pd.DataFrame(
        {
            "date": dates,
            "value1": np.cumsum(np.random.randn(100)),
            "value2": np.cumsum(np.random.randn(100)) * 2,
            "value3": np.sin(np.linspace(0, 4 * np.pi, 100)) * 10,
        }
    )
    ts_file = test_dir / "time_series.csv"
    ts_data.to_csv(ts_file, index=False)
    print(f"  ✓ Created time series: {ts_file}")

    # 3. Mixed categorical/numeric
    mixed_data = pd.DataFrame(
        {
            "category": np.random.choice(["A", "B", "C", "D"], 200),
            "group": np.random.choice(["Control", "Treatment"], 200),
            "age": np.random.randint(20, 80, 200),
            "score": np.random.normal(50, 15, 200),
            "measurement": np.random.exponential(2, 200),
        }
    )
    mixed_file = test_dir / "mixed_data.csv"
    mixed_data.to_csv(mixed_file, index=False)
    print(f"  ✓ Created mixed data: {mixed_file}")

    # 4. Genomic intervals
    genomic_data = pd.DataFrame(
        {
            "chr": [f"chr{i}" for i in np.random.randint(1, 23, 50)],
            "start": np.random.randint(1000, 100000, 50),
            "end": np.random.randint(100000, 200000, 50),
            "score": np.random.uniform(0, 100, 50),
            "strand": np.random.choice(["+", "-"], 50),
        }
    )
    genomic_file = test_dir / "genomic_intervals.csv"
    genomic_data.to_csv(genomic_file, index=False)
    print(f"  ✓ Created genomic intervals: {genomic_file}")

    return test_dir, {
        "expression": expr_file,
        "timeseries": ts_file,
        "mixed": mixed_file,
        "genomic": genomic_file,
    }


def test_analyze(file_path):
    """Test file analysis endpoint"""
    print(f"\n🔍 Analyzing: {file_path.name}")

    response = requests.post(
        f"{API_BASE_URL}/api/quick-visualizer/analyze",
        json={"file_path": str(file_path)},
    )

    if response.status_code != 200:
        print(f"  ❌ Failed: {response.status_code} - {response.text}")
        return None

    result = response.json()

    if not result.get("success"):
        print(f"  ❌ Analysis failed: {result.get('error')}")
        return None

    print(f"  ✓ Data type: {result['data_type']}")
    print(f"  ✓ Dimensions: {result['n_rows']} rows × {result['n_columns']} columns")

    if result.get("suggestions"):
        print(f"  ✓ Suggestions ({len(result['suggestions'])}):")
        for i, sug in enumerate(result["suggestions"][:3], 1):
            print(f"    {i}. {sug['title']}: {sug['description']}")

    return result


def test_visualize(file_path, chart_type, options=None):
    """Test visualization generation endpoint"""
    print(f"\n📊 Generating {chart_type} for: {file_path.name}")

    response = requests.post(
        f"{API_BASE_URL}/api/quick-visualizer/visualize",
        json={
            "file_path": str(file_path),
            "chart_type": chart_type,
            "options": options or {},
        },
    )

    if response.status_code != 200:
        print(f"  ❌ Failed: {response.status_code} - {response.text}")
        return False

    result = response.json()

    if not result.get("success"):
        print(f"  ❌ Visualization failed: {result.get('error')}")
        return False

    # Check image data
    image_data = result.get("image", "")
    if image_data.startswith("data:image/png;base64,"):
        image_size = len(image_data)
        print(f"  ✓ Generated PNG image ({image_size} bytes)")
        return True
    else:
        print(f"  ❌ Invalid image format")
        return False


def test_list_chart_types():
    """Test chart types listing endpoint"""
    print(f"\n📋 Listing chart types")

    response = requests.get(f"{API_BASE_URL}/api/quick-visualizer/chart-types")

    if response.status_code != 200:
        print(f"  ❌ Failed: {response.status_code}")
        return

    chart_types = response.json()
    print(f"  ✓ Found {len(chart_types)} chart types:")
    for ct in chart_types:
        print(f"    - {ct['type']}: {ct['name']}")


def test_list_data_types():
    """Test data types listing endpoint"""
    print(f"\n📋 Listing data types")

    response = requests.get(f"{API_BASE_URL}/api/quick-visualizer/data-types")

    if response.status_code != 200:
        print(f"  ❌ Failed: {response.status_code}")
        return

    result = response.json()
    data_types = result.get("data_types", [])
    print(f"  ✓ Found {len(data_types)} data types:")
    for dt in data_types:
        print(f"    - {dt['type']}: {dt['name']}")


def run_all_tests():
    """Run all tests"""
    print("=" * 60)
    print("🧪 Quick Visualization API Tests")
    print("=" * 60)

    # Test metadata endpoints
    test_list_chart_types()
    test_list_data_types()

    # Create test data
    test_dir, test_files = create_test_data()

    try:
        # Test analysis for all files
        results = {}
        for name, file_path in test_files.items():
            results[name] = test_analyze(file_path)

        # Test visualizations based on suggestions
        print("\n" + "=" * 60)
        print("📊 Testing Visualizations")
        print("=" * 60)

        # Expression matrix - heatmap
        if results.get("expression"):
            test_visualize(test_files["expression"], "heatmap")
            test_visualize(test_files["expression"], "correlation")

        # Time series - line chart
        if results.get("timeseries"):
            test_visualize(test_files["timeseries"], "line")
            test_visualize(
                test_files["timeseries"],
                "scatter",
                {"x_column": "value1", "y_column": "value2"},
            )

        # Mixed data - various charts
        if results.get("mixed"):
            test_visualize(test_files["mixed"], "histogram", {"column": "age"})
            test_visualize(test_files["mixed"], "bar", {"column": "category"})
            test_visualize(
                test_files["mixed"], "scatter", {"x_column": "age", "y_column": "score"}
            )
            test_visualize(
                test_files["mixed"], "boxplot", {"columns": ["age", "score"]}
            )

        # Genomic intervals
        if results.get("genomic"):
            test_visualize(test_files["genomic"], "histogram", {"column": "score"})

        print("\n" + "=" * 60)
        print("✅ All tests completed!")
        print("=" * 60)

    finally:
        # Cleanup
        import shutil

        shutil.rmtree(test_dir)
        print(f"\n🧹 Cleaned up test data: {test_dir}")


if __name__ == "__main__":
    run_all_tests()
