#!/usr/bin/env Rscript
#
# Bioconverter R Wrapper Package
# ===============================
#
# Complete R wrapper for bioconverter Python package.
# Provides convenient R interface to bioconverter functionality.
#
# Usage:
#   source("bioconverter_wrapper.R")
#   bc_init()
#   bc_convert("input.vcf", "output.csv")
#

library(reticulate)

# Package-level variables
.bioconverter_env <- new.env(parent = emptyenv())
.bioconverter_env$module <- NULL
.bioconverter_env$initialized <- FALSE


#' Initialize bioconverter module
#'
#' @param force Force re-initialization
#' @return TRUE if successful
#' @export
bc_init <- function(force = FALSE) {
  
  if (.bioconverter_env$initialized && !force) {
    return(invisible(TRUE))
  }
  
  tryCatch({
    .bioconverter_env$module <- import("bioconverter")
    .bioconverter_env$initialized <- TRUE
    
    version <- .bioconverter_env$module$`__version__`
    message(sprintf("Bioconverter initialized (version %s)", version))
    
    return(invisible(TRUE))
    
  }, error = function(e) {
    warning("Failed to initialize bioconverter. Install with: py_install('bioconverter')")
    return(invisible(FALSE))
  })
}


#' Get bioconverter module
#'
#' @return bioconverter Python module
.get_module <- function() {
  if (!.bioconverter_env$initialized) {
    bc_init()
  }
  return(.bioconverter_env$module)
}


#' Convert file format
#'
#' @param input_file Input file path
#' @param output_file Output file path
#' @param from_format Source format (optional, auto-detected)
#' @param to_format Target format (optional, inferred from extension)
#' @param validate Validate input file before conversion
#' @param ... Additional arguments passed to converter
#' @return TRUE if successful, FALSE otherwise
#' @export
#'
#' @examples
#' bc_convert("input.vcf", "output.csv")
#' bc_convert("reads.fastq", "reads.fasta", validate = TRUE)
bc_convert <- function(input_file, output_file, 
                       from_format = NULL, 
                       to_format = NULL,
                       validate = FALSE,
                       ...) {
  
  module <- .get_module()
  converter <- module$Converter()
  
  # Build arguments
  args <- list(
    input_file = input_file,
    output_file = output_file,
    validate = validate
  )
  
  if (!is.null(from_format)) args$from_format <- from_format
  if (!is.null(to_format)) args$to_format <- to_format
  
  # Add additional arguments
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    args <- c(args, extra_args)
  }
  
  tryCatch({
    do.call(converter$convert, args)
    return(invisible(TRUE))
    
  }, error = function(e) {
    warning(sprintf("Conversion failed: %s", e$message))
    return(invisible(FALSE))
  })
}


#' List available converters
#'
#' @param format_type Filter by format type (optional)
#' @return Data frame with available converters
#' @export
#'
#' @examples
#' bc_list_converters()
#' bc_list_converters(format_type = "vcf")
bc_list_converters <- function(format_type = NULL) {
  
  module <- .get_module()
  
  tryCatch({
    if (is.null(format_type)) {
      converters <- module$list_converters()
    } else {
      converters <- module$list_converters(format_type = format_type)
    }
    
    # Convert to data frame for easier viewing in R
    if (length(converters) > 0) {
      df <- do.call(rbind, lapply(converters, function(x) {
        data.frame(
          from = x$from,
          to = x$to,
          stringsAsFactors = FALSE
        )
      }))
      return(df)
    } else {
      return(data.frame(from = character(), to = character()))
    }
    
  }, error = function(e) {
    warning(sprintf("Failed to list converters: %s", e$message))
    return(NULL)
  })
}


#' Detect file format
#'
#' @param file_path File path or content to detect
#' @return Detected format string
#' @export
#'
#' @examples
#' bc_detect_format("input.vcf")
bc_detect_format <- function(file_path) {
  
  module <- .get_module()
  
  tryCatch({
    format <- module$detect_format(file_path)
    return(format)
    
  }, error = function(e) {
    warning(sprintf("Format detection failed: %s", e$message))
    return(NULL)
  })
}


#' Validate file format
#'
#' @param file_path File path to validate
#' @param format_type Format type (optional, auto-detected)
#' @return List with is_valid (logical) and errors (character vector)
#' @export
#'
#' @examples
#' result <- bc_validate("input.vcf")
#' if (result$is_valid) {
#'   message("File is valid")
#' }
bc_validate <- function(file_path, format_type = NULL) {
  
  module <- .get_module()
  
  # Auto-detect format if not provided
  if (is.null(format_type)) {
    format_type <- bc_detect_format(file_path)
    if (is.null(format_type)) {
      return(list(is_valid = FALSE, errors = "Could not detect format"))
    }
  }
  
  tryCatch({
    result <- module$validate_file(file_path, format_type = format_type)
    
    return(list(
      is_valid = result[[1]],
      errors = if (length(result[[2]]) > 0) result[[2]] else character(0)
    ))
    
  }, error = function(e) {
    return(list(
      is_valid = FALSE,
      errors = e$message
    ))
  })
}


#' Get bioconverter version
#'
#' @return Version string
#' @export
bc_version <- function() {
  module <- .get_module()
  return(module$`__version__`)
}


#' Batch convert files
#'
#' @param input_files Vector of input file paths
#' @param output_dir Output directory
#' @param to_format Target format
#' @param verbose Print progress
#' @return Data frame with conversion results
#' @export
#'
#' @examples
#' files <- list.files("data", pattern = "\\.vcf$", full.names = TRUE)
#' bc_batch_convert(files, "output", to_format = "csv")
bc_batch_convert <- function(input_files, output_dir, 
                              to_format = "csv",
                              verbose = TRUE) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize results
  results <- data.frame(
    input_file = character(),
    output_file = character(),
    success = logical(),
    error = character(),
    stringsAsFactors = FALSE
  )
  
  # Progress bar
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(input_files), style = 3)
  }
  
  # Convert each file
  for (i in seq_along(input_files)) {
    input_file <- input_files[i]
    base_name <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(output_dir, sprintf("%s.%s", base_name, to_format))
    
    success <- bc_convert(input_file, output_file)
    
    results <- rbind(results, data.frame(
      input_file = input_file,
      output_file = output_file,
      success = success,
      error = ifelse(success, "", "Conversion failed"),
      stringsAsFactors = FALSE
    ))
    
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) {
    close(pb)
    message(sprintf("\nCompleted: %d/%d successful", 
                   sum(results$success), 
                   nrow(results)))
  }
  
  return(results)
}


#' Convert and load into data frame
#'
#' @param file_path Input file path
#' @param to_format Intermediate format (default: "csv")
#' @param ... Additional arguments passed to read function
#' @return Data frame with file contents
#' @export
#'
#' @examples
#' data <- bc_read("variants.vcf")
#' head(data)
bc_read <- function(file_path, to_format = "csv", ...) {
  
  # Create temporary file
  temp_file <- tempfile(fileext = sprintf(".%s", to_format))
  
  # Convert
  success <- bc_convert(file_path, temp_file)
  
  if (!success) {
    stop("Conversion failed")
  }
  
  # Read based on format
  data <- tryCatch({
    if (to_format == "csv") {
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fread(temp_file, ...)
      } else {
        read.csv(temp_file, ...)
      }
    } else if (to_format == "tsv") {
      read.delim(temp_file, ...)
    } else {
      stop(sprintf("Unsupported format for reading: %s", to_format))
    }
  }, error = function(e) {
    stop(sprintf("Failed to read converted file: %s", e$message))
  })
  
  # Clean up
  unlink(temp_file)
  
  return(data)
}


# Print package information when sourced
if (!exists(".bioconverter_loaded")) {
  message("
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║                    Bioconverter R Wrapper                            ║
║                                                                      ║
║  Available functions:                                                ║
║    bc_init()              - Initialize bioconverter                  ║
║    bc_convert()           - Convert file format                      ║
║    bc_list_converters()   - List available converters                ║
║    bc_detect_format()     - Detect file format                       ║
║    bc_validate()          - Validate file format                     ║
║    bc_version()           - Get version                              ║
║    bc_batch_convert()     - Batch convert files                      ║
║    bc_read()              - Convert and read into data frame         ║
║                                                                      ║
║  Example:                                                            ║
║    bc_init()                                                         ║
║    bc_convert('input.vcf', 'output.csv')                             ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
  ")
  
  .bioconverter_loaded <- TRUE
}
