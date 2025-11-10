#!/usr/bin/env Rscript
#
# Data Frame Integration Example
# ===============================
#
# Demonstrates converting bioinformatics files and loading into R data frames.
#
# Usage:
#     Rscript dataframe_integration.R <input_vcf>
#
# Example:
#     Rscript dataframe_integration.R sample.vcf
#

# Load required libraries
if (!requireNamespace("reticulate", quietly = TRUE)) {
  message("Error: reticulate package is required")
  message("Install it with: install.packages('reticulate')")
  quit(status = 1)
}

library(reticulate)

# Try to load data.table for faster reading
has_datatable <- requireNamespace("data.table", quietly = TRUE)
if (has_datatable) {
  library(data.table)
}


#' Convert VCF to CSV and load into R data frame
#'
#' @param vcf_file Path to VCF file
#' @param keep_csv Keep the intermediate CSV file (default: FALSE)
#' @return Data frame with VCF data
convert_and_load <- function(vcf_file, keep_csv = FALSE) {
  
  # Check if input file exists
  if (!file.exists(vcf_file)) {
    stop(sprintf("Error: Input file '%s' not found", vcf_file))
  }
  
  # Import bioconverter
  tryCatch({
    bioconverter <- import("bioconverter")
  }, error = function(e) {
    stop("bioconverter is not installed. Install with: py_install('bioconverter')")
  })
  
  message(sprintf("Processing %s...\n", vcf_file))
  
  # Create temporary CSV file
  temp_csv <- tempfile(fileext = ".csv")
  
  # Convert VCF to CSV
  message("Step 1: Converting VCF to CSV...")
  
  tryCatch({
    converter <- bioconverter$Converter()
    converter$convert(vcf_file, temp_csv)
    message("  ✓ Conversion successful\n")
  }, error = function(e) {
    stop(sprintf("Conversion failed: %s", e$message))
  })
  
  # Load into R
  message("Step 2: Loading data into R...")
  
  data <- tryCatch({
    if (has_datatable) {
      fread(temp_csv, showProgress = FALSE)
    } else {
      read.csv(temp_csv, stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    stop(sprintf("Failed to load CSV: %s", e$message))
  })
  
  message(sprintf("  ✓ Loaded %d rows and %d columns\n", nrow(data), ncol(data)))
  
  # Keep or remove temporary CSV
  if (keep_csv) {
    final_csv <- sub("\\.vcf$", ".csv", vcf_file)
    file.copy(temp_csv, final_csv, overwrite = TRUE)
    message(sprintf("CSV saved to: %s\n", final_csv))
  } else {
    unlink(temp_csv)
  }
  
  return(data)
}


#' Analyze VCF data in R
#'
#' @param data VCF data frame
#' @return Summary statistics
analyze_vcf_data <- function(data) {
  
  message(strrep("=", 60))
  message("Data Analysis")
  message(strrep("=", 60), "\n")
  
  # Basic statistics
  message("Dataset dimensions:")
  message(sprintf("  Rows: %d", nrow(data)))
  message(sprintf("  Columns: %d", ncol(data)))
  message(sprintf("\nColumn names: %s", paste(names(data), collapse = ", ")))
  
  # Show first few rows
  message("\nFirst 3 rows:")
  print(head(data, 3))
  
  # Summary statistics for numeric columns
  numeric_cols <- sapply(data, is.numeric)
  if (any(numeric_cols)) {
    message("\nSummary of numeric columns:")
    print(summary(data[, numeric_cols, drop = FALSE]))
  }
  
  # Quality filtering example (if QUAL column exists)
  if ("QUAL" %in% names(data)) {
    message("\nQuality filtering example:")
    high_quality <- data[data$QUAL > 30, ]
    message(sprintf("  Variants with QUAL > 30: %d (%.1f%%)", 
                   nrow(high_quality), 
                   100 * nrow(high_quality) / nrow(data)))
  }
  
  message(sprintf("\n%s\n", strrep("=", 60)))
  
  return(invisible(data))
}


#' Export analysis results
#'
#' @param data Data frame to export
#' @param output_file Output file path
export_results <- function(data, output_file = "analysis_results.csv") {
  
  message(sprintf("Exporting results to %s...", output_file))
  
  tryCatch({
    write.csv(data, output_file, row.names = FALSE)
    
    size_kb <- file.info(output_file)$size / 1024
    message(sprintf("  ✓ Exported %d rows (%.2f KB)\n", nrow(data), size_kb))
    
  }, error = function(e) {
    message(sprintf("  ✗ Export failed: %s\n", e$message))
  })
}


# Main execution
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check arguments
  if (length(args) < 1) {
    cat("Usage: Rscript dataframe_integration.R <input_vcf>\n\n")
    cat("Example:\n")
    cat("  Rscript dataframe_integration.R sample.vcf\n")
    cat("  Rscript dataframe_integration.R variants.vcf\n")
    quit(status = 1)
  }
  
  vcf_file <- args[1]
  
  # Convert and load
  tryCatch({
    data <- convert_and_load(vcf_file, keep_csv = FALSE)
    
    # Analyze
    analyze_vcf_data(data)
    
    # Optional: Filter and export
    # Example: Export high-quality variants
    if ("QUAL" %in% names(data)) {
      filtered_data <- data[data$QUAL > 30, ]
      export_results(filtered_data, "high_quality_variants.csv")
    }
    
  }, error = function(e) {
    message(sprintf("Error: %s", e$message))
    quit(status = 1)
  })
  
  message("Processing complete!")
  quit(status = 0)
}

# If running interactively, provide usage example
if (interactive()) {
  message("
Example usage in interactive R session:

  library(reticulate)
  source('dataframe_integration.R')
  
  # Convert and load VCF data
  vcf_data <- convert_and_load('sample.vcf')
  
  # Analyze the data
  analyze_vcf_data(vcf_data)
  
  # Filter data
  high_qual <- vcf_data[vcf_data$QUAL > 30, ]
  
  # Export results
  export_results(high_qual, 'filtered_variants.csv')
  
  # Or from command line:
  # Rscript dataframe_integration.R sample.vcf
  ")
}
