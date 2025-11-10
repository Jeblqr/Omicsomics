# Bioconverter Usage Guide

## Overview

**Bioconverter** is a Python package for converting between various bioinformatics file formats. It supports genomics, transcriptomics, proteomics, and metabolomics data formats with a simple, unified API.

**PyPI Package**: https://pypi.org/project/bioconverter/0.1.0/

---

## Table of Contents

- [Installation](#installation)
- [Python Usage](#python-usage)
  - [Command Line Interface (CLI)](#python-cli)
  - [Interactive Python](#interactive-python)
  - [Python Scripts](#python-scripts)
- [R Usage](#r-usage)
  - [Command Line Interface (CLI)](#r-cli)
  - [Interactive R](#interactive-r)
  - [R Scripts](#r-scripts)
- [Supported Formats](#supported-formats)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

---

## Installation

### Python Installation

**Requirements**: Python 3.8+

#### Install from PyPI

```bash
# Install the latest version
pip install bioconverter

# Install specific version
pip install bioconverter==0.1.0

# Upgrade to latest version
pip install --upgrade bioconverter
```

#### Install in Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv bioconv-env

# Activate virtual environment
# On Linux/Mac:
source bioconv-env/bin/activate
# On Windows:
bioconv-env\Scripts\activate

# Install bioconverter
pip install bioconverter
```

#### Install with Additional Dependencies

```bash
# Install with all optional dependencies
pip install bioconverter[all]

# Install with specific extras
pip install bioconverter[genomics]
pip install bioconverter[proteomics]
```

#### Verify Installation

```bash
# Check version
bioconverter --version

# Show help
bioconverter --help

# List available converters
bioconverter list
```

### R Installation

**Requirements**: R 4.0+ and Python 3.8+

#### Install via reticulate

```r
# Install reticulate if not already installed
install.packages("reticulate")

# Load reticulate
library(reticulate)

# Install bioconverter in Python
py_install("bioconverter")

# Verify installation
bioconverter <- import("bioconverter")
bioconverter$`__version__`
```

#### Using Conda Environment (Recommended)

```r
# Create conda environment
conda_create("bioconv-r", python_version = "3.10")

# Install bioconverter in the environment
use_condaenv("bioconv-r")
py_install("bioconverter", envname = "bioconv-r")

# Use the environment
library(reticulate)
use_condaenv("bioconv-r")
bioconverter <- import("bioconverter")
```

---

## Python Usage

### Python CLI

#### Basic Conversion

```bash
# Convert VCF to CSV
bioconverter convert input.vcf output.csv

# Convert with explicit format specification
bioconverter convert input.vcf output.csv --from-format vcf --to-format csv

# Convert FASTQ to FASTA
bioconverter convert reads.fastq reads.fasta
```

#### List Available Converters

```bash
# List all converters
bioconverter list

# List converters for specific format
bioconverter list --format vcf

# Show detailed information
bioconverter list --verbose
```

#### Get Format Information

```bash
# Show info about a format
bioconverter info vcf

# Show supported conversions
bioconverter info vcf --conversions
```

#### Batch Conversion

```bash
# Convert multiple files
bioconverter convert *.vcf --output-dir converted/

# Convert with pattern matching
bioconverter convert input/*.fastq --to-format fasta --output-dir fasta_output/
```

#### Advanced Options

```bash
# Validate input file before conversion
bioconverter convert input.vcf output.csv --validate

# Compression
bioconverter convert input.vcf output.csv.gz --compress

# Set quality threshold for FASTQ
bioconverter convert input.fastq output.fasta --quality-threshold 30

# Include metadata in output
bioconverter convert input.vcf output.json --include-metadata
```

### Interactive Python

#### Basic Usage

```python
from bioconverter import Converter

# Initialize converter
converter = Converter()

# Convert VCF to CSV
converter.convert(
    input_file="input.vcf",
    output_file="output.csv"
)

# Convert with explicit formats
converter.convert(
    input_file="input.vcf",
    output_file="output.csv",
    from_format="vcf",
    to_format="csv"
)
```

#### List Available Converters

```python
from bioconverter import list_converters

# Get all converters
converters = list_converters()
for conv in converters:
    print(f"{conv['from']} -> {conv['to']}")

# Get converters for specific format
vcf_converters = list_converters(format_type="vcf")
print(vcf_converters)
```

#### Format Detection

```python
from bioconverter import detect_format

# Auto-detect format
file_format = detect_format("input.vcf")
print(f"Detected format: {file_format}")

# Detect from content
with open("input.vcf", "r") as f:
    content = f.read(1000)  # Read first 1KB
    file_format = detect_format(content=content)
```

#### Validation

```python
from bioconverter import validate_file

# Validate file format
is_valid, errors = validate_file("input.vcf", format_type="vcf")

if is_valid:
    print("File is valid")
else:
    print("Validation errors:")
    for error in errors:
        print(f"  - {error}")
```

### Python Scripts

#### Example 1: Simple Conversion Script

```python
#!/usr/bin/env python3
"""
Simple VCF to CSV conversion script
"""
from bioconverter import Converter
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: convert_vcf.py <input.vcf> <output.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Create converter
    converter = Converter()

    # Convert
    try:
        converter.convert(input_file, output_file)
        print(f"Successfully converted {input_file} to {output_file}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
```

#### Example 2: Batch Conversion with Progress

```python
#!/usr/bin/env python3
"""
Batch convert all VCF files in a directory to CSV
"""
from bioconverter import Converter
from pathlib import Path
from tqdm import tqdm

def batch_convert(input_dir, output_dir, from_format="vcf", to_format="csv"):
    """Convert all files in input_dir to output_dir"""

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Find all input files
    files = list(input_path.glob(f"*.{from_format}"))

    if not files:
        print(f"No .{from_format} files found in {input_dir}")
        return

    # Convert with progress bar
    converter = Converter()

    for input_file in tqdm(files, desc="Converting"):
        output_file = output_path / f"{input_file.stem}.{to_format}"

        try:
            converter.convert(str(input_file), str(output_file))
        except Exception as e:
            print(f"\nError converting {input_file}: {e}")

if __name__ == "__main__":
    batch_convert("vcf_files/", "csv_output/")
```

#### Example 3: Conversion with Validation

```python
#!/usr/bin/env python3
"""
Convert with validation and error handling
"""
from bioconverter import Converter, validate_file, detect_format
import sys

def safe_convert(input_file, output_file):
    """Convert with validation and error handling"""

    # Detect format
    print(f"Detecting format of {input_file}...")
    detected_format = detect_format(input_file)
    print(f"Detected format: {detected_format}")

    # Validate
    print("Validating file...")
    is_valid, errors = validate_file(input_file, format_type=detected_format)

    if not is_valid:
        print("Validation failed:")
        for error in errors:
            print(f"  - {error}")
        return False

    print("File is valid")

    # Convert
    print(f"Converting to {output_file}...")
    converter = Converter()

    try:
        converter.convert(
            input_file,
            output_file,
            validate=True,
            include_metadata=True
        )
        print("Conversion successful!")
        return True
    except Exception as e:
        print(f"Conversion failed: {e}", file=sys.stderr)
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: safe_convert.py <input> <output>")
        sys.exit(1)

    success = safe_convert(sys.argv[1], sys.argv[2])
    sys.exit(0 if success else 1)
```

#### Example 4: Using Context Manager

```python
from bioconverter import Converter

# Using context manager for automatic cleanup
with Converter() as converter:
    # Convert multiple files
    files = ["file1.vcf", "file2.vcf", "file3.vcf"]

    for input_file in files:
        output_file = input_file.replace(".vcf", ".csv")
        converter.convert(input_file, output_file)
        print(f"Converted {input_file}")
```

#### Example 5: Programmatic Format Conversion

```python
from bioconverter import Converter

# Initialize converter
converter = Converter()

# Genomics conversions
converter.convert("variants.vcf", "variants.csv")
converter.convert("reads.fastq", "reads.fasta")
converter.convert("genes.gff", "genes.bed")

# Proteomics conversions
converter.convert("peptides.mzML", "peptides.csv")
converter.convert("proteins.fasta", "proteins.csv")

# Metabolomics conversions
converter.convert("metabolites.mzXML", "metabolites.csv")

# Get conversion statistics
stats = converter.get_stats()
print(f"Converted {stats['total_files']} files")
print(f"Success rate: {stats['success_rate']:.1%}")
```

---

## R Usage

### R CLI

Since bioconverter is a Python package, R users can call it via command line:

```r
# In R, use system() to call CLI
system("bioconverter convert input.vcf output.csv")

# With output capture
result <- system("bioconverter convert input.vcf output.csv", intern = TRUE)
print(result)

# Check if conversion succeeded
exit_code <- system("bioconverter convert input.vcf output.csv")
if (exit_code == 0) {
  message("Conversion successful")
} else {
  stop("Conversion failed")
}
```

### Interactive R

#### Setup and Basic Usage

```r
# Load reticulate
library(reticulate)

# Import bioconverter
bioconverter <- import("bioconverter")

# Create converter
converter <- bioconverter$Converter()

# Convert file
converter$convert(
  input_file = "input.vcf",
  output_file = "output.csv"
)

# Convert with explicit formats
converter$convert(
  input_file = "input.vcf",
  output_file = "output.csv",
  from_format = "vcf",
  to_format = "csv"
)
```

#### List Available Converters

```r
library(reticulate)
bioconverter <- import("bioconverter")

# Get all converters
converters <- bioconverter$list_converters()

# Print converter information
for (i in seq_along(converters)) {
  conv <- converters[[i]]
  cat(sprintf("%s -> %s\n", conv$from, conv$to))
}

# Get converters for specific format
vcf_converters <- bioconverter$list_converters(format_type = "vcf")
print(vcf_converters)
```

#### Format Detection

```r
library(reticulate)
bioconverter <- import("bioconverter")

# Detect format
file_format <- bioconverter$detect_format("input.vcf")
cat("Detected format:", file_format, "\n")

# Detect from content
content <- readLines("input.vcf", n = 100)
content_str <- paste(content, collapse = "\n")
file_format <- bioconverter$detect_format(content = content_str)
```

#### Validation

```r
library(reticulate)
bioconverter <- import("bioconverter")

# Validate file
validation <- bioconverter$validate_file("input.vcf", format_type = "vcf")
is_valid <- validation[[1]]
errors <- validation[[2]]

if (is_valid) {
  message("File is valid")
} else {
  message("Validation errors:")
  for (error in errors) {
    message("  - ", error)
  }
}
```

### R Scripts

#### Example 1: Simple Conversion Function

```r
#!/usr/bin/env Rscript
#
# Simple VCF to CSV conversion function
#

library(reticulate)

convert_vcf_to_csv <- function(input_file, output_file) {
  # Import bioconverter
  bioconverter <- import("bioconverter")

  # Create converter
  converter <- bioconverter$Converter()

  # Convert
  tryCatch({
    converter$convert(input_file, output_file)
    message(sprintf("Successfully converted %s to %s", input_file, output_file))
    return(TRUE)
  }, error = function(e) {
    message(sprintf("Error: %s", e$message))
    return(FALSE)
  })
}

# Usage
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    cat("Usage: convert_vcf.R <input.vcf> <output.csv>\n")
    quit(status = 1)
  }

  success <- convert_vcf_to_csv(args[1], args[2])
  quit(status = ifelse(success, 0, 1))
}
```

#### Example 2: Batch Conversion in R

```r
#!/usr/bin/env Rscript
#
# Batch convert all VCF files to CSV
#

library(reticulate)

batch_convert <- function(input_dir, output_dir, from_format = "vcf", to_format = "csv") {
  # Import bioconverter
  bioconverter <- import("bioconverter")
  converter <- bioconverter$Converter()

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Find all input files
  pattern <- sprintf("\\.%s$", from_format)
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    message(sprintf("No .%s files found in %s", from_format, input_dir))
    return(invisible())
  }

  # Convert files with progress
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)

  for (i in seq_along(files)) {
    input_file <- files[i]
    base_name <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(output_dir, sprintf("%s.%s", base_name, to_format))

    tryCatch({
      converter$convert(input_file, output_file)
    }, error = function(e) {
      message(sprintf("\nError converting %s: %s", input_file, e$message))
    })

    setTxtProgressBar(pb, i)
  }

  close(pb)
  message(sprintf("\nConverted %d files", length(files)))
}

# Usage
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    cat("Usage: batch_convert.R <input_dir> <output_dir>\n")
    quit(status = 1)
  }

  batch_convert(args[1], args[2])
}
```

#### Example 3: Integration with R Data Frames

```r
#!/usr/bin/env Rscript
#
# Convert and load into R data frame
#

library(reticulate)
library(data.table)

convert_and_load <- function(vcf_file) {
  # Import bioconverter
  bioconverter <- import("bioconverter")
  converter <- bioconverter$Converter()

  # Create temporary CSV file
  temp_csv <- tempfile(fileext = ".csv")

  # Convert VCF to CSV
  tryCatch({
    converter$convert(vcf_file, temp_csv)
    message("Conversion successful")
  }, error = function(e) {
    stop(sprintf("Conversion failed: %s", e$message))
  })

  # Load into R
  data <- fread(temp_csv)

  # Clean up
  unlink(temp_csv)

  return(data)
}

# Usage example
vcf_data <- convert_and_load("variants.vcf")
head(vcf_data)
summary(vcf_data)

# Perform R analysis
filtered_data <- vcf_data[vcf_data$QUAL > 30, ]
nrow(filtered_data)
```

#### Example 4: R Package Wrapper

```r
#' Bioconverter R Wrapper
#'
#' Provides convenient R interface to bioconverter
#'
#' @export

library(reticulate)

# Initialize module
.bioconverter <- NULL

#' Initialize bioconverter
#'
#' @export
init_bioconverter <- function() {
  if (is.null(.bioconverter)) {
    .bioconverter <<- import("bioconverter")
  }
  invisible(.bioconverter)
}

#' Convert file format
#'
#' @param input_file Input file path
#' @param output_file Output file path
#' @param from_format Source format (optional)
#' @param to_format Target format (optional)
#' @param validate Validate input file
#' @return TRUE if successful
#' @export
bc_convert <- function(input_file, output_file, from_format = NULL,
                       to_format = NULL, validate = FALSE) {
  init_bioconverter()
  converter <- .bioconverter$Converter()

  args <- list(
    input_file = input_file,
    output_file = output_file,
    validate = validate
  )

  if (!is.null(from_format)) args$from_format <- from_format
  if (!is.null(to_format)) args$to_format <- to_format

  tryCatch({
    do.call(converter$convert, args)
    return(TRUE)
  }, error = function(e) {
    warning(sprintf("Conversion failed: %s", e$message))
    return(FALSE)
  })
}

#' List available converters
#'
#' @param format_type Filter by format type
#' @return List of converters
#' @export
bc_list_converters <- function(format_type = NULL) {
  init_bioconverter()

  if (is.null(format_type)) {
    .bioconverter$list_converters()
  } else {
    .bioconverter$list_converters(format_type = format_type)
  }
}

#' Detect file format
#'
#' @param file_path File path
#' @return Detected format
#' @export
bc_detect_format <- function(file_path) {
  init_bioconverter()
  .bioconverter$detect_format(file_path)
}

#' Validate file
#'
#' @param file_path File path
#' @param format_type Format type
#' @return List with is_valid and errors
#' @export
bc_validate <- function(file_path, format_type = NULL) {
  init_bioconverter()

  if (is.null(format_type)) {
    format_type <- bc_detect_format(file_path)
  }

  result <- .bioconverter$validate_file(file_path, format_type = format_type)

  list(
    is_valid = result[[1]],
    errors = result[[2]]
  )
}

# Usage examples
# bc_convert("input.vcf", "output.csv")
# converters <- bc_list_converters()
# format <- bc_detect_format("file.vcf")
# validation <- bc_validate("file.vcf")
```

---

## Supported Formats

### Genomics

| Format  | Extension   | Description                  |
| ------- | ----------- | ---------------------------- |
| VCF     | .vcf        | Variant Call Format          |
| FASTQ   | .fastq, .fq | Sequence with quality scores |
| FASTA   | .fasta, .fa | Sequence format              |
| SAM     | .sam        | Sequence Alignment/Map       |
| BAM     | .bam        | Binary SAM                   |
| BED     | .bed        | Browser Extensible Data      |
| GFF/GTF | .gff, .gtf  | Gene annotation              |
| MAF     | .maf        | Mutation Annotation Format   |

### Transcriptomics

| Format       | Extension  | Description                    |
| ------------ | ---------- | ------------------------------ |
| Count Matrix | .txt, .csv | Gene expression counts         |
| TPM          | .txt, .csv | Transcripts Per Million        |
| FPKM         | .txt, .csv | Fragments Per Kilobase Million |

### Proteomics

| Format | Extension  | Description              |
| ------ | ---------- | ------------------------ |
| mzML   | .mzml      | Mass spectrometry data   |
| mzXML  | .mzxml     | Mass spectrometry XML    |
| MGF    | .mgf       | Mascot Generic Format    |
| PSM    | .txt, .csv | Peptide-Spectrum Matches |

### Metabolomics

| Format | Extension | Description            |
| ------ | --------- | ---------------------- |
| mzTab  | .mztab    | Metabolomics data      |
| mzData | .mzdata   | Mass spectrometry data |

### Universal

| Format  | Extension | Description                |
| ------- | --------- | -------------------------- |
| CSV     | .csv      | Comma-separated values     |
| TSV     | .tsv      | Tab-separated values       |
| JSON    | .json     | JavaScript Object Notation |
| Excel   | .xlsx     | Microsoft Excel            |
| Parquet | .parquet  | Apache Parquet             |

---

## Examples

### Example 1: VCF to CSV

**Python**:

```python
from bioconverter import Converter

converter = Converter()
converter.convert("variants.vcf", "variants.csv")
```

**R**:

```r
library(reticulate)
bioconverter <- import("bioconverter")
converter <- bioconverter$Converter()
converter$convert("variants.vcf", "variants.csv")
```

**CLI**:

```bash
bioconverter convert variants.vcf variants.csv
```

### Example 2: FASTQ to FASTA

**Python**:

```python
from bioconverter import Converter

converter = Converter()
converter.convert(
    "reads.fastq",
    "reads.fasta",
    quality_threshold=30
)
```

**R**:

```r
library(reticulate)
bioconverter <- import("bioconverter")
converter <- bioconverter$Converter()
converter$convert(
  input_file = "reads.fastq",
  output_file = "reads.fasta",
  quality_threshold = 30L
)
```

**CLI**:

```bash
bioconverter convert reads.fastq reads.fasta --quality-threshold 30
```

### Example 3: Batch Conversion

**Python**:

```python
from bioconverter import Converter
from pathlib import Path

converter = Converter()
vcf_files = Path("data").glob("*.vcf")

for vcf_file in vcf_files:
    csv_file = vcf_file.with_suffix(".csv")
    converter.convert(str(vcf_file), str(csv_file))
```

**R**:

```r
library(reticulate)
bioconverter <- import("bioconverter")
converter <- bioconverter$Converter()

vcf_files <- list.files("data", pattern = "\\.vcf$", full.names = TRUE)

for (vcf_file in vcf_files) {
  csv_file <- sub("\\.vcf$", ".csv", vcf_file)
  converter$convert(vcf_file, csv_file)
}
```

**CLI**:

```bash
bioconverter convert data/*.vcf --output-dir csv_output/
```

---

## Troubleshooting

### Common Issues

#### 1. Import Error in Python

```python
# Error: ModuleNotFoundError: No module named 'bioconverter'

# Solution: Install bioconverter
pip install bioconverter

# Or check if you're using the correct Python environment
which python
pip list | grep bioconverter
```

#### 2. Import Error in R

```r
# Error: Module 'bioconverter' not found

# Solution: Install in Python used by reticulate
library(reticulate)
py_install("bioconverter")

# Or specify Python version
use_python("/usr/bin/python3")
py_install("bioconverter")
```

#### 3. Conversion Failed

```python
# Error: Unsupported conversion

# Solution: Check supported conversions
from bioconverter import list_converters
converters = list_converters()
print(converters)

# Or use explicit format specification
converter.convert(
    "input.txt",
    "output.csv",
    from_format="custom",
    to_format="csv"
)
```

#### 4. File Not Found

```bash
# Error: FileNotFoundError: [Errno 2] No such file or directory

# Solution: Use absolute path or check current directory
pwd
ls -la input.vcf

# Or use absolute path
bioconverter convert /full/path/to/input.vcf output.csv
```

#### 5. Permission Denied

```bash
# Error: PermissionError: [Errno 13] Permission denied

# Solution: Check file permissions
ls -la input.vcf

# Or change permissions
chmod 644 input.vcf

# Or run with sudo (not recommended)
sudo bioconverter convert input.vcf output.csv
```

### Getting Help

```bash
# Command line help
bioconverter --help
bioconverter convert --help
bioconverter list --help

# Python help
python -c "from bioconverter import Converter; help(Converter)"

# R help
library(reticulate)
bioconverter <- import("bioconverter")
help(bioconverter$Converter)
```

### Reporting Issues

If you encounter bugs or have feature requests:

1. **Check existing issues**: Search the issue tracker
2. **Provide details**: Include version, OS, error messages
3. **Minimal example**: Provide code to reproduce the issue
4. **Submit issue**: Create a new issue on the repository

```bash
# Get version information
bioconverter --version
python --version
pip show bioconverter
```

---

## Additional Resources

- **PyPI Package**: https://pypi.org/project/bioconverter/0.1.0/
- **Documentation**: Check package documentation with `help(bioconverter)`
- **Examples**: See `examples/` directory in source repository
- **Issue Tracker**: Report bugs and request features

---

_Last Updated: 2025-11-10_
