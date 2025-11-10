# Bioconverter Quick Reference

> Quick reference card for bioconverter - bioinformatics format conversion tool

## 📦 Installation

```bash
# Python
pip install bioconverter

# R
library(reticulate)
py_install("bioconverter")
```

## 🔧 Command Line (CLI)

```bash
# Basic conversion
bioconverter convert input.vcf output.csv

# Explicit formats
bioconverter convert input.vcf output.csv --from-format vcf --to-format csv

# List converters
bioconverter list
bioconverter list --format vcf

# Show format info
bioconverter info vcf

# Batch conversion
bioconverter convert *.vcf --output-dir csv_output/

# With validation
bioconverter convert input.vcf output.csv --validate

# Compression
bioconverter convert input.vcf output.csv.gz --compress
```

## 🐍 Python

### Basic Usage

```python
from bioconverter import Converter

converter = Converter()
converter.convert("input.vcf", "output.csv")
```

### List Converters

```python
from bioconverter import list_converters

converters = list_converters()
vcf_converters = list_converters(format_type="vcf")
```

### Format Detection

```python
from bioconverter import detect_format

format_type = detect_format("input.vcf")
```

### Validation

```python
from bioconverter import validate_file

is_valid, errors = validate_file("input.vcf", format_type="vcf")
if not is_valid:
    for error in errors:
        print(f"Error: {error}")
```

### Batch Conversion

```python
from bioconverter import Converter
from pathlib import Path

converter = Converter()
for vcf_file in Path("data").glob("*.vcf"):
    csv_file = vcf_file.with_suffix(".csv")
    converter.convert(str(vcf_file), str(csv_file))
```

### Context Manager

```python
from bioconverter import Converter

with Converter() as converter:
    converter.convert("file1.vcf", "file1.csv")
    converter.convert("file2.vcf", "file2.csv")
```

## 📊 R

### Basic Usage

```r
library(reticulate)
bioconverter <- import("bioconverter")

converter <- bioconverter$Converter()
converter$convert("input.vcf", "output.csv")
```

### List Converters

```r
converters <- bioconverter$list_converters()
vcf_converters <- bioconverter$list_converters(format_type = "vcf")
```

### Format Detection

```r
format_type <- bioconverter$detect_format("input.vcf")
```

### Validation

```r
validation <- bioconverter$validate_file("input.vcf", format_type = "vcf")
is_valid <- validation[[1]]
errors <- validation[[2]]
```

### Batch Conversion

```r
converter <- bioconverter$Converter()
files <- list.files("data", pattern = "\\.vcf$", full.names = TRUE)

for (file in files) {
  output <- sub("\\.vcf$", ".csv", file)
  converter$convert(file, output)
}
```

### Using R Wrapper

```r
source("bioconverter_wrapper.R")
bc_init()

# Convert
bc_convert("input.vcf", "output.csv")

# List converters
bc_list_converters()

# Detect format
bc_detect_format("input.vcf")

# Validate
bc_validate("input.vcf")

# Read into data frame
data <- bc_read("input.vcf")
```

## 📁 Supported Formats

### Genomics

VCF, FASTQ, FASTA, SAM, BAM, BED, GFF, GTF, MAF

### Transcriptomics

Count Matrix, TPM, FPKM

### Proteomics

mzML, mzXML, MGF, PSM

### Metabolomics

mzTab, mzData

### Universal

CSV, TSV, JSON, Excel, Parquet

## 🔍 Common Conversions

```bash
# Genomics
bioconverter convert variants.vcf variants.csv
bioconverter convert reads.fastq reads.fasta
bioconverter convert genes.gff genes.bed

# Proteomics
bioconverter convert peptides.mzML peptides.csv

# Universal
bioconverter convert data.csv data.json
bioconverter convert data.csv data.xlsx
```

## 🐛 Troubleshooting

### Module Not Found

```bash
# Python
pip install bioconverter
pip list | grep bioconverter

# R
library(reticulate)
py_install("bioconverter")
```

### File Not Found

```bash
# Use absolute path
bioconverter convert /full/path/to/input.vcf output.csv

# Or check current directory
pwd
ls -la input.vcf
```

### Check Version

```bash
bioconverter --version
```

```python
import bioconverter
print(bioconverter.__version__)
```

```r
bioconverter <- import("bioconverter")
bioconverter$`__version__`
```

## 📚 Resources

- **PyPI**: https://pypi.org/project/bioconverter/
- **Documentation**: `docs/BIOCONVERTER_USAGE.md`
- **Examples**: `examples/bioconverter/`
- **Tests**: `scripts/test_bioconverter.py`

## 💡 Pro Tips

1. **Use validation** for important conversions
2. **Batch convert** for efficiency
3. **Check formats** with `list_converters()`
4. **Auto-detect** formats when possible
5. **Use context managers** for cleanup
6. **Read directly** into data frames (R)
7. **Pipeline** conversions for multi-step transforms

## 🎯 Example Workflows

### Python: VCF Analysis Pipeline

```python
from bioconverter import Converter

converter = Converter()

# 1. Convert VCF to CSV
converter.convert("variants.vcf", "variants.csv")

# 2. Convert CSV to JSON for web
converter.convert("variants.csv", "variants.json")

# 3. Convert to Excel for sharing
converter.convert("variants.csv", "variants.xlsx")
```

### R: Load and Analyze VCF

```r
library(reticulate)
source("bioconverter_wrapper.R")
bc_init()

# Convert and load in one step
vcf_data <- bc_read("variants.vcf")

# Analyze
high_qual <- vcf_data[vcf_data$QUAL > 30, ]

# Export results
write.csv(high_qual, "filtered_variants.csv")
```

### CLI: Batch Genomics Pipeline

```bash
# Convert all VCF files
bioconverter convert *.vcf --output-dir csv/

# Convert all FASTQ to FASTA
bioconverter convert *.fastq --to-format fasta --output-dir fasta/

# Convert GFF to BED
bioconverter convert genes.gff genes.bed
```

---

_For detailed documentation, see `docs/BIOCONVERTER_USAGE.md`_
