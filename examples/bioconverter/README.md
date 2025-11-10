# Bioconverter Examples

This directory contains practical examples for using the bioconverter package.

## 📁 Contents

### Python Examples

- `simple_convert.py` - Basic file conversion
- `batch_convert.py` - Batch conversion with progress tracking
- `validate_and_convert.py` - Conversion with validation
- `advanced_usage.py` - Advanced features demonstration

### R Examples

- `simple_convert.R` - Basic file conversion in R
- `batch_convert.R` - Batch conversion in R
- `dataframe_integration.R` - Load converted data into R data frames
- `bioconverter_wrapper.R` - Complete R package wrapper

### Test Data

- Sample input files for testing conversions

## 🚀 Quick Start

### Python

```bash
# Run simple conversion
python simple_convert.py ../test_data/sample_genomics.vcf output.csv

# Run batch conversion
python batch_convert.py ../test_data/ output_dir/

# Run with validation
python validate_and_convert.py ../test_data/sample_genomics.vcf output.csv
```

### R

```bash
# Run simple conversion
Rscript simple_convert.R ../test_data/sample_genomics.vcf output.csv

# Run batch conversion
Rscript batch_convert.R ../test_data/ output_dir/
```

## 📚 Documentation

See `docs/BIOCONVERTER_USAGE.md` for comprehensive usage guide.

## 🧪 Testing

All examples are designed to work with the test data in `test_data/` directory.
