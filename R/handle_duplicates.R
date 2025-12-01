#' Remove Duplicate Gene Symbols
#'
#' Handle duplicate gene symbols using different strategies
#'
#' @param df Data frame with SYMBOL column
#' @param method Character. Method for handling duplicates: "random", "average", "highest", "first"
#' @param symbol_col Character. Name of symbol column (default: "SYMBOL")
#' @param count_cols Character vector. Names of count columns (if NULL, auto-detect)
#' @return Data frame with unique gene symbols
#' @export
#' @examples
#' # Random selection (fastest)
#' unique_counts <- remove_duplicates(counts, method = "random")
#' 
#' # Average expression (conservative)
#' unique_counts <- remove_duplicates(counts, method = "average")
#' 
#' # Keep highest expressing isoform
#' unique_counts <- remove_duplicates(counts, method = "highest")

library(dplyr)

remove_duplicates <- function(df, 
                             method = c("random", "average", "highest", "first"),
                             symbol_col = "SYMBOL",
                             count_cols = NULL) {
  
  method <- match.arg(method)
  
  if (!symbol_col %in% colnames(df)) {
    stop(paste("Column", symbol_col, "not found"))
  }
  
  # Identify count columns if not specified
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  # Check for duplicates
  n_before <- nrow(df)
  n_unique <- length(unique(df[[symbol_col]]))
  n_duplicates <- n_before - n_unique
  
  cat("Total genes:", n_before, "\n")
  cat("Unique symbols:", n_unique, "\n")
  cat("Duplicate entries:", n_duplicates, "\n")
  
  if (n_duplicates == 0) {
    cat("No duplicates found. Returning original data.\n")
    return(df)
  }
  
  cat("Using method:", method, "\n")
  
  # Apply selected method
  if (method == "random") {
    # Randomly select one entry per symbol (fastest)
    df_unique <- df[!duplicated(df[[symbol_col]]), ]
    cat("Randomly selected one entry per duplicate symbol\n")
    
  } else if (method == "first") {
    # Keep first occurrence
    df_unique <- df[!duplicated(df[[symbol_col]]), ]
    cat("Kept first occurrence of each duplicate symbol\n")
    
  } else if (method == "average") {
    # Average expression across duplicates
    # Separate count and non-count columns
    non_count_cols <- setdiff(colnames(df), count_cols)
    
    # For count columns: calculate mean
    df_counts <- df[, c(symbol_col, count_cols)]
    df_counts_avg <- aggregate(. ~ get(symbol_col), data = df_counts, FUN = mean)
    colnames(df_counts_avg)[1] <- symbol_col
    
    # For non-count columns: take first occurrence
    df_meta <- df[!duplicated(df[[symbol_col]]), non_count_cols, drop = FALSE]
    
    # Merge back
    df_unique <- merge(df_meta, df_counts_avg, by = symbol_col)
    cat("Averaged expression across duplicate symbols\n")
    
  } else if (method == "highest") {
    # Keep isoform with highest mean expression
    df$mean_expression <- rowMeans(df[, count_cols, drop = FALSE])
    df <- df[order(df[[symbol_col]], -df$mean_expression), ]
    df_unique <- df[!duplicated(df[[symbol_col]]), ]
    df_unique$mean_expression <- NULL
    cat("Kept highest expressing isoform per duplicate symbol\n")
  }
  
  cat("\nFinal gene count:", nrow(df_unique), "\n")
  
  # Set symbols as rownames
  rownames(df_unique) <- df_unique[[symbol_col]]
  
  return(df_unique)
}


#' Identify Duplicate Genes
#'
#' Find which gene symbols appear multiple times
#'
#' @param df Data frame with SYMBOL column
#' @param symbol_col Character. Name of symbol column
#' @return Data frame showing duplicate symbols and their counts
#' @export

find_duplicates <- function(df, symbol_col = "SYMBOL") {
  
  if (!symbol_col %in% colnames(df)) {
    stop(paste("Column", symbol_col, "not found"))
  }
  
  symbol_counts <- table(df[[symbol_col]])
  duplicated_symbols <- symbol_counts[symbol_counts > 1]
  
  if (length(duplicated_symbols) == 0) {
    cat("No duplicate symbols found.\n")
    return(NULL)
  }
  
  # Create summary dataframe
  dup_df <- data.frame(
    Symbol = names(duplicated_symbols),
    Count = as.numeric(duplicated_symbols)
  )
  dup_df <- dup_df[order(-dup_df$Count), ]
  
  cat("Found", nrow(dup_df), "symbols with duplicates\n")
  cat("Total duplicate entries:", sum(dup_df$Count - 1), "\n")
  cat("\nTop duplicated symbols:\n")
  print(head(dup_df, 10))
  
  return(dup_df)
}


#' Compare Duplicate Handling Methods
#'
#' Show statistics for different duplicate handling approaches
#'
#' @param df Data frame with SYMBOL column and count data
#' @param symbol_col Character. Name of symbol column
#' @return List with comparison statistics
#' @export

compare_duplicate_methods <- function(df, symbol_col = "SYMBOL") {
  
  methods <- c("random", "average", "highest")
  results <- list()
  
  cat("Comparing duplicate handling methods...\n\n")
  
  for (method in methods) {
    cat("Method:", method, "\n")
    df_processed <- remove_duplicates(df, method = method, symbol_col = symbol_col)
    
    results[[method]] <- list(
      n_genes = nrow(df_processed),
      method = method
    )
    cat("\n")
  }
  
  return(results)
}
