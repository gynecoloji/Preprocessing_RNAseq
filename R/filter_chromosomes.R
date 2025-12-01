#' Filter Chromosomes
#'
#' Remove genes from unwanted chromosomes (Y, MT, unplaced contigs, etc.)
#'
#' @param df Data frame with Chr column
#' @param keep_pattern Character. Regex pattern for chromosomes to keep (default: autosomes + X)
#' @param exclude_pattern Character. Regex pattern for chromosomes to exclude (optional)
#' @param chr_col Character. Name of chromosome column (default: "Chr")
#' @return Filtered data frame
#' @export
#' @examples
#' # Keep only standard chromosomes (1-22, X)
#' filtered <- filter_chromosomes(counts, keep_pattern = "^[0-9X]+")
#' 
#' # Exclude Y, MT, and unplaced contigs
#' filtered <- filter_chromosomes(counts, exclude_pattern = "^(Y|MT|KI|GL)")

library(stringr)

filter_chromosomes <- function(df, 
                              keep_pattern = "^[0-9X]+",
                              exclude_pattern = NULL,
                              chr_col = "Chr") {
  
  if (!chr_col %in% colnames(df)) {
    stop(paste("Column", chr_col, "not found in data frame"))
  }
  
  n_before <- nrow(df)
  cat("Genes before chromosome filtering:", n_before, "\n")
  
  # Apply keep pattern
  if (!is.null(keep_pattern)) {
    df <- df[grep(keep_pattern, df[[chr_col]]), ]
    cat("After keeping pattern '", keep_pattern, "':", nrow(df), "genes\n", sep = "")
  }
  
  # Apply exclude pattern
  if (!is.null(exclude_pattern)) {
    exclude_idx <- grep(exclude_pattern, df[[chr_col]])
    if (length(exclude_idx) > 0) {
      df <- df[-exclude_idx, ]
      cat("After excluding pattern '", exclude_pattern, "':", nrow(df), "genes\n", sep = "")
    }
  }
  
  # Check chromosome distribution
  # Split multi-chromosome annotations (some genes span multiple chr)
  chr_split <- str_split(df[[chr_col]], ";", n = 2, simplify = TRUE)
  chr_table <- table(chr_split[, 1])
  
  cat("\nChromosome distribution:\n")
  print(chr_table)
  
  n_removed <- n_before - nrow(df)
  cat("\nRemoved", n_removed, "genes (", 
      round(n_removed/n_before * 100, 2), "%)\n", sep = "")
  
  return(df)
}


#' Check Chromosome Distribution
#'
#' Summarize the distribution of genes across chromosomes
#'
#' @param df Data frame with Chr column
#' @param chr_col Character. Name of chromosome column
#' @return Table of chromosome counts
#' @export

check_chromosome_distribution <- function(df, chr_col = "Chr") {
  
  if (!chr_col %in% colnames(df)) {
    stop(paste("Column", chr_col, "not found in data frame"))
  }
  
  # Handle multi-chromosome annotations
  chr_split <- str_split(df[[chr_col]], ";", n = 2, simplify = TRUE)
  chr_table <- table(chr_split[, 1])
  
  # Sort by chromosome number
  chr_names <- names(chr_table)
  numeric_chr <- suppressWarnings(as.numeric(chr_names))
  
  # Separate numeric and non-numeric chromosomes
  is_numeric <- !is.na(numeric_chr)
  
  if (any(is_numeric)) {
    sorted_numeric <- names(sort(chr_table[is_numeric]))
    sorted_other <- names(sort(chr_table[!is_numeric]))
    chr_table <- chr_table[c(sorted_numeric, sorted_other)]
  }
  
  return(chr_table)
}
