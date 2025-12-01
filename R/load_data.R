#' Load FeatureCounts Output
#'
#' Reads count data from FeatureCounts output file and returns a cleaned count matrix
#'
#' @param file_path Character. Path to the FeatureCounts output file
#' @param skip_comment Character. Comment character to skip (default: "#")
#' @return Data frame with gene metadata and count columns
#' @export
#' @examples
#' counts <- load_featurecounts("data/featureCount.txt")

load_featurecounts <- function(file_path, skip_comment = "#") {
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Read the data
  cat("Loading FeatureCounts data from:", file_path, "\n")
  
  df <- read.table(
    file_path,
    sep = "\t",
    header = TRUE,
    comment.char = skip_comment,
    row.names = 1,
    check.names = FALSE
  )
  
  cat("Loaded", nrow(df), "genes and", ncol(df), "columns\n")
  
  # Show first few rows and columns
  cat("\nFirst 5 rows and columns:\n")
  print(df[1:min(5, nrow(df)), 1:min(5, ncol(df))])
  
  return(df)
}


#' Get Count Matrix from FeatureCounts Output
#'
#' Extract only the count columns from FeatureCounts data
#'
#' @param df Data frame from load_featurecounts()
#' @param metadata_cols Character vector. Column names to exclude (default: standard metadata)
#' @return Data frame with only count data
#' @export

get_count_matrix <- function(df, 
                             metadata_cols = c("Chr", "Start", "End", "Strand", "Length")) {
  
  # Remove metadata columns
  count_cols <- setdiff(colnames(df), metadata_cols)
  count_matrix <- df[, count_cols, drop = FALSE]
  
  cat("Extracted count matrix with", nrow(count_matrix), "genes and", 
      ncol(count_matrix), "samples\n")
  
  return(count_matrix)
}


#' Get Metadata from FeatureCounts Output
#'
#' Extract gene metadata (Chr, Start, End, etc.) from FeatureCounts data
#'
#' @param df Data frame from load_featurecounts()
#' @param metadata_cols Character vector. Column names to extract
#' @return Data frame with metadata only
#' @export

get_metadata <- function(df, 
                        metadata_cols = c("Chr", "Start", "End", "Strand", "Length")) {
  
  # Keep only metadata columns that exist
  available_cols <- intersect(metadata_cols, colnames(df))
  
  if (length(available_cols) == 0) {
    warning("No metadata columns found")
    return(NULL)
  }
  
  metadata <- df[, available_cols, drop = FALSE]
  metadata$ENSEMBL_ID <- rownames(df)
  
  return(metadata)
}
