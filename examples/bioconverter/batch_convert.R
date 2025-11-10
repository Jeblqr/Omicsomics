#!/usr/bin/env Rscript
#
# Batch Bioconverter Example for R
# =================================
#
# Converts all files of a specific format in a directory.
#
# Usage:
#     Rscript batch_convert.R <input_dir> <output_dir> [from_format] [to_format]
#
# Example:
#     Rscript batch_convert.R vcf_files/ csv_output/
#     Rscript batch_convert.R data/ converted/ vcf csv
#

# Load required library
if (!requireNamespace("reticulate", quietly = TRUE)) {
  message("Error: reticulate package is required")
  message("Install it with: install.packages('reticulate')")
  quit(status = 1)
}

library(reticulate)


#' Batch convert files
#'
#' @param input_dir Input directory path
#' @param output_dir Output directory path
#' @param from_format Source file format (default: "vcf")
#' @param to_format Target file format (default: "csv")
#' @param verbose Print progress information (default: TRUE)
#' @return List with conversion statistics
batch_convert <- function(input_dir, output_dir, 
                          from_format = "vcf", 
                          to_format = "csv",
                          verbose = TRUE) {
  
  # Check if input directory exists
  if (!dir.exists(input_dir)) {
    message(sprintf("Error: Input directory '%s' not found", input_dir))
    return(invisible(NULL))
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Import bioconverter
  tryCatch({
    bioconverter <- import("bioconverter")
  }, error = function(e) {
    message("Error: bioconverter is not installed in Python")
    message("\nInstall it with:")
    message("  library(reticulate)")
    message("  py_install('bioconverter')")
    return(invisible(NULL))
  })
  
  # Create converter
  converter <- bioconverter$Converter()
  
  # Find all input files
  pattern <- sprintf("\\.%s$", from_format)
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    message(sprintf("No .%s files found in %s", from_format, input_dir))
    return(invisible(NULL))
  }
  
  if (verbose) {
    message(sprintf("Found %d .%s files", length(files), from_format))
    message(sprintf("Converting to .%s format...\n", to_format))
  }
  
  # Initialize statistics
  stats <- list(
    total = length(files),
    success = 0,
    errors = 0,
    error_details = list()
  )
  
  # Create progress bar
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  }
  
  # Convert files
  for (i in seq_along(files)) {
    input_file <- files[i]
    base_name <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(output_dir, sprintf("%s.%s", base_name, to_format))
    
    tryCatch({
      converter$convert(input_file, output_file)
      stats$success <- stats$success + 1
      
      if (verbose && !exists("pb")) {
        message(sprintf("  ✓ %s -> %s", basename(input_file), basename(output_file)))
      }
      
    }, error = function(e) {
      stats$errors <<- stats$errors + 1
      stats$error_details[[length(stats$error_details) + 1]] <<- list(
        file = basename(input_file),
        error = e$message
      )
      
      if (verbose) {
        message(sprintf("\n  ✗ Error converting %s: %s", basename(input_file), e$message))
      }
    })
    
    if (verbose && exists("pb")) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose && exists("pb")) {
    close(pb)
  }
  
  # Print summary
  if (verbose) {
    message(sprintf("\n%s", strrep("=", 60)))
    message("Conversion Summary:")
    message(sprintf("  Success: %d/%d", stats$success, stats$total))
    message(sprintf("  Errors:  %d/%d", stats$errors, stats$total))
    message(sprintf("  Output directory: %s", normalizePath(output_dir)))
    message(strrep("=", 60))
    
    # Show error details if any
    if (stats$errors > 0 && length(stats$error_details) > 0) {
      message("\nError details:")
      for (err in stats$error_details) {
        message(sprintf("  - %s: %s", err$file, substr(err$error, 1, 50)))
      }
    }
  }
  
  return(invisible(stats))
}


# Main execution
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check arguments
  if (length(args) < 2) {
    cat("Usage: Rscript batch_convert.R <input_dir> <output_dir> [from_format] [to_format]\n\n")
    cat("Example:\n")
    cat("  Rscript batch_convert.R vcf_files/ csv_output/\n")
    cat("  Rscript batch_convert.R data/ converted/ vcf csv\n")
    cat("  Rscript batch_convert.R fastq/ fasta/ fastq fasta\n")
    quit(status = 1)
  }
  
  input_dir <- args[1]
  output_dir <- args[2]
  from_format <- ifelse(length(args) >= 3, args[3], "vcf")
  to_format <- ifelse(length(args) >= 4, args[4], "csv")
  
  # Run batch conversion
  stats <- batch_convert(
    input_dir,
    output_dir,
    from_format = from_format,
    to_format = to_format,
    verbose = TRUE
  )
  
  # Exit with appropriate status
  exit_code <- ifelse(is.null(stats) || stats$errors > 0, 1, 0)
  quit(status = exit_code)
}

# If running interactively, provide usage example
if (interactive()) {
  message("
Example usage in interactive R session:

  library(reticulate)
  source('batch_convert.R')
  
  # Batch convert VCF to CSV
  batch_convert('vcf_files/', 'csv_output/')
  
  # Convert FASTQ to FASTA
  batch_convert('fastq/', 'fasta/', from_format='fastq', to_format='fasta')
  
  # Or from command line:
  # Rscript batch_convert.R vcf_files/ csv_output/ vcf csv
  ")
}
