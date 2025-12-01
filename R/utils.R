#' Rename Sample Columns
#'
#' Rename count matrix columns to more readable sample names
#'
#' @param df Data frame with count data
#' @param old_names Character vector. Current column names
#' @param new_names Character vector. New column names
#' @param auto_detect Logical. Auto-detect count columns (default: TRUE)
#' @return Data frame with renamed columns
#' @export
#' @examples
#' # Manual renaming
#' renamed <- rename_samples(counts, 
#'                          old_names = c("sample1.bam", "sample2.bam"),
#'                          new_names = c("Control_1", "Treatment_1"))

rename_samples <- function(df, 
                          old_names = NULL, 
                          new_names = NULL,
                          auto_detect = TRUE) {
  
  if (is.null(old_names) || is.null(new_names)) {
    stop("Both old_names and new_names must be provided")
  }
  
  if (length(old_names) != length(new_names)) {
    stop("old_names and new_names must have the same length")
  }
  
  # Check if old names exist
  missing_names <- setdiff(old_names, colnames(df))
  if (length(missing_names) > 0) {
    warning(paste("These columns not found:", paste(missing_names, collapse = ", ")))
  }
  
  # Rename
  for (i in seq_along(old_names)) {
    if (old_names[i] %in% colnames(df)) {
      colnames(df)[colnames(df) == old_names[i]] <- new_names[i]
    }
  }
  
  cat("Renamed", length(old_names), "columns\n")
  cat("New column names:\n")
  print(colnames(df))
  
  return(df)
}


#' Save Filtered Count Matrix
#'
#' Save count matrix with metadata about filtering steps
#'
#' @param df Data frame to save
#' @param output_file Character. Output file path
#' @param compression Logical. Compress output (default: FALSE)
#' @param include_metadata Logical. Include processing metadata (default: TRUE)
#' @return Invisible NULL
#' @export

save_count_matrix <- function(df, 
                              output_file,
                              compression = FALSE,
                              include_metadata = TRUE) {
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save main file
  write.csv(df, output_file, row.names = TRUE)
  cat("Saved count matrix to:", output_file, "\n")
  cat("Dimensions:", nrow(df), "genes x", ncol(df), "columns\n")
  
  # Save metadata file
  if (include_metadata) {
    metadata_file <- sub("\\.csv$", "_metadata.txt", output_file)
    
    sink(metadata_file)
    cat("RNA-seq Count Matrix Metadata\n")
    cat("==============================\n\n")
    cat("File:", output_file, "\n")
    cat("Date:", as.character(Sys.time()), "\n")
    cat("Dimensions:", nrow(df), "genes x", ncol(df), "columns\n\n")
    
    cat("Column names:\n")
    cat(paste(colnames(df), collapse = ", "), "\n\n")
    
    cat("Summary statistics:\n")
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
    
    if (length(count_cols) > 0) {
      count_matrix <- df[, count_cols, drop = FALSE]
      cat("Total counts per sample:\n")
      print(colSums(count_matrix))
      cat("\nMean counts per sample:\n")
      print(colMeans(count_matrix))
    }
    
    sink()
    
    cat("Saved metadata to:", metadata_file, "\n")
  }
  
  invisible(NULL)
}


#' Load Preprocessed Count Matrix
#'
#' Load a previously saved count matrix
#'
#' @param file_path Character. Path to CSV file
#' @return Data frame with count data
#' @export

load_count_matrix <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  df <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  cat("Loaded count matrix from:", file_path, "\n")
  cat("Dimensions:", nrow(df), "genes x", ncol(df), "columns\n")
  
  return(df)
}


#' Get Sample Information
#'
#' Extract sample names and basic info from count matrix
#'
#' @param df Data frame with count data
#' @param count_cols Character vector. Names of count columns (if NULL, auto-detect)
#' @return Data frame with sample information
#' @export

get_sample_info <- function(df, count_cols = NULL) {
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  count_matrix <- df[, count_cols, drop = FALSE]
  
  sample_df <- data.frame(
    Sample = colnames(count_matrix),
    TotalCounts = colSums(count_matrix),
    MeanCounts = colMeans(count_matrix),
    MedianCounts = apply(count_matrix, 2, median),
    GenesDetected = colSums(count_matrix > 0)
  )
  
  rownames(sample_df) <- NULL
  
  return(sample_df)
}


#' Print Processing Summary
#'
#' Print summary of preprocessing steps
#'
#' @param original_df Original data frame (before processing)
#' @param final_df Final data frame (after processing)
#' @return Invisible NULL
#' @export

print_processing_summary <- function(original_df, final_df) {
  
  cat("\n")
  cat("=" , rep("=", 50), "\n", sep = "")
  cat("PREPROCESSING SUMMARY\n")
  cat("=" , rep("=", 50), "\n", sep = "")
  
  cat("\nOriginal data:\n")
  cat("  Genes:", nrow(original_df), "\n")
  
  metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
  orig_count_cols <- setdiff(colnames(original_df), metadata_cols)
  cat("  Samples:", length(orig_count_cols), "\n")
  
  cat("\nFinal data:\n")
  cat("  Genes:", nrow(final_df), "\n")
  
  final_count_cols <- setdiff(colnames(final_df), metadata_cols)
  cat("  Samples:", length(final_count_cols), "\n")
  
  cat("\nGenes removed:", nrow(original_df) - nrow(final_df), 
      "(", round((nrow(original_df) - nrow(final_df))/nrow(original_df) * 100, 1), "%)\n", sep = "")
  
  cat("\nFinal sample names:\n")
  cat("  ", paste(final_count_cols, collapse = ", "), "\n", sep = "")
  
  cat("=" , rep("=", 50), "\n", sep = "")
  cat("\n")
  
  invisible(NULL)
}
