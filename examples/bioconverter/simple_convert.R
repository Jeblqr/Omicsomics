#!/usr/bin/env Rscript
#
# Simple Bioconverter Example for R
# ==================================
#
# Demonstrates basic file format conversion using bioconverter in R.
#
# Usage:
#     Rscript simple_convert.R <input_file> <output_file>
#
# Example:
#     Rscript simple_convert.R input.vcf output.csv
#

# Load required library
if (!requireNamespace("reticulate", quietly = TRUE)) {
  message("Error: reticulate package is required")
  message("Install it with: install.packages('reticulate')")
  quit(status = 1)
}

library(reticulate)


#' Simple file conversion function
#'
#' @param input_file Path to input file
#' @param output_file Path to output file
#' @return TRUE if successful, FALSE otherwise
simple_convert <- function(input_file, output_file) {
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    message(sprintf("Error: Input file '%s' not found", input_file))
    return(FALSE)
  }
  
  # Import bioconverter
  tryCatch({
    bioconverter <- import("bioconverter")
  }, error = function(e) {
    message("Error: bioconverter is not installed in Python")
    message("\nInstall it with:")
    message("  library(reticulate)")
    message("  py_install('bioconverter')")
    return(FALSE)
  })
  
  # Create converter
  converter <- bioconverter$Converter()
  
  # Perform conversion
  message(sprintf("Converting %s to %s...", input_file, output_file))
  
  tryCatch({
    converter$convert(input_file, output_file)
    
    message(sprintf("✓ Successfully converted to %s", output_file))
    
    # Show output file size
    if (file.exists(output_file)) {
      size_kb <- file.info(output_file)$size / 1024
      message(sprintf("  Output size: %.2f KB", size_kb))
    }
    
    return(TRUE)
    
  }, error = function(e) {
    message(sprintf("Error during conversion: %s", e$message))
    return(FALSE)
  })
}


# Main execution
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check arguments
  if (length(args) != 2) {
    cat("Usage: Rscript simple_convert.R <input_file> <output_file>\n\n")
    cat("Example:\n")
    cat("  Rscript simple_convert.R input.vcf output.csv\n")
    cat("  Rscript simple_convert.R reads.fastq reads.fasta\n")
    quit(status = 1)
  }
  
  input_file <- args[1]
  output_file <- args[2]
  
  # Run conversion
  success <- simple_convert(input_file, output_file)
  
  # Exit with appropriate status
  quit(status = ifelse(success, 0, 1))
}

# If running interactively, provide usage example
if (interactive()) {
  message("
Example usage in interactive R session:

  library(reticulate)
  source('simple_convert.R')
  
  # Convert file
  simple_convert('input.vcf', 'output.csv')
  
  # Or from command line:
  # Rscript simple_convert.R input.vcf output.csv
  ")
}
