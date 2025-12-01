#' Filter Low Expression Genes
#'
#' Remove genes with low expression across samples
#'
#' @param df Data frame with count data
#' @param min_mean Numeric. Minimum mean expression across all samples (default: 10)
#' @param min_count Numeric. Minimum count per sample (alternative method)
#' @param min_samples Integer. Minimum number of samples meeting min_count (alternative method)
#' @param count_cols Character vector. Names of count columns (if NULL, auto-detect)
#' @return Filtered data frame
#' @export
#' @examples
#' # Filter by mean expression
#' filtered <- filter_low_expression(counts, min_mean = 10)
#' 
#' # Filter by counts in minimum samples (e.g., CPM > 5 in at least 3 samples)
#' filtered <- filter_low_expression(counts, min_count = 5, min_samples = 3)

filter_low_expression <- function(df, 
                                  min_mean = NULL,
                                  min_count = NULL,
                                  min_samples = NULL,
                                  count_cols = NULL) {
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  count_matrix <- df[, count_cols, drop = FALSE]
  
  n_before <- nrow(df)
  cat("Genes before filtering:", n_before, "\n")
  
  # Method 1: Filter by mean expression
  if (!is.null(min_mean)) {
    gene_means <- rowMeans(count_matrix)
    keep <- gene_means > min_mean
    
    cat("Filtering by mean expression > ", min_mean, "\n", sep = "")
    cat("Genes passing filter:", sum(keep), "\n")
    
    df <- df[keep, ]
  }
  
  # Method 2: Filter by minimum count in minimum samples
  if (!is.null(min_count) && !is.null(min_samples)) {
    samples_passing <- rowSums(count_matrix > min_count)
    keep <- samples_passing >= min_samples
    
    cat("Filtering by count > ", min_count, " in at least ", 
        min_samples, " samples\n", sep = "")
    cat("Genes passing filter:", sum(keep), "\n")
    
    df <- df[keep, ]
  }
  
  n_removed <- n_before - nrow(df)
  cat("\nRemoved:", n_removed, "genes (", 
      round(n_removed/n_before * 100, 2), "%)\n", sep = "")
  cat("Remaining:", nrow(df), "genes\n")
  
  return(df)
}


#' Filter by Gene Type
#'
#' Keep only specific gene types (e.g., protein-coding)
#'
#' @param df Data frame with GENETYPE column
#' @param gene_type Character or vector. Gene type(s) to keep (default: "protein-coding")
#' @param genetype_col Character. Name of gene type column (default: "GENETYPE")
#' @param add_genetype Logical. Add GENETYPE column if missing (default: TRUE)
#' @param organism Character. Organism database for annotation if needed
#' @return Filtered data frame
#' @export
#' @examples
#' # Keep only protein-coding genes
#' protein_coding <- filter_by_gene_type(counts, gene_type = "protein-coding")
#' 
#' # Keep multiple gene types
#' filtered <- filter_by_gene_type(counts, 
#'                                 gene_type = c("protein-coding", "lncRNA"))

filter_by_gene_type <- function(df, 
                                gene_type = "protein-coding",
                                genetype_col = "GENETYPE",
                                add_genetype = TRUE,
                                organism = "org.Hs.eg.db") {
  
  # Add GENETYPE column if missing
  if (!genetype_col %in% colnames(df) && add_genetype) {
    cat("GENETYPE column not found. Adding annotation...\n")
    source("R/annotate_genes.R")
    df <- annotate_gene_types(df, organism = organism)
  }
  
  if (!genetype_col %in% colnames(df)) {
    stop(paste("Column", genetype_col, "not found and could not be added"))
  }
  
  n_before <- nrow(df)
  cat("Genes before filtering:", n_before, "\n")
  
  # Show gene type distribution
  cat("\nGene type distribution:\n")
  print(table(df[[genetype_col]], useNA = "ifany"))
  
  # Filter
  keep <- df[[genetype_col]] %in% gene_type
  df_filtered <- df[keep, ]
  
  cat("\nKeeping gene type(s):", paste(gene_type, collapse = ", "), "\n")
  cat("Genes passing filter:", nrow(df_filtered), "\n")
  
  n_removed <- n_before - nrow(df_filtered)
  cat("Removed:", n_removed, "genes (", 
      round(n_removed/n_before * 100, 2), "%)\n", sep = "")
  
  return(df_filtered)
}


#' Get Expression Statistics
#'
#' Calculate summary statistics for expression levels
#'
#' @param df Data frame with count data
#' @param count_cols Character vector. Names of count columns (if NULL, auto-detect)
#' @return Data frame with statistics per gene
#' @export

get_expression_stats <- function(df, count_cols = NULL) {
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  count_matrix <- df[, count_cols, drop = FALSE]
  
  stats_df <- data.frame(
    gene = rownames(count_matrix),
    mean = rowMeans(count_matrix),
    median = apply(count_matrix, 1, median),
    sd = apply(count_matrix, 1, sd),
    min = apply(count_matrix, 1, min),
    max = apply(count_matrix, 1, max),
    n_zero = rowSums(count_matrix == 0),
    n_nonzero = rowSums(count_matrix > 0)
  )
  
  stats_df <- stats_df[order(-stats_df$mean), ]
  
  return(stats_df)
}


#' Plot Expression Distribution
#'
#' Visualize distribution of expression levels
#'
#' @param df Data frame with count data
#' @param count_cols Character vector. Names of count columns
#' @param log_transform Logical. Log transform counts (default: TRUE)
#' @return ggplot object
#' @export

plot_expression_distribution <- function(df, 
                                        count_cols = NULL,
                                        log_transform = TRUE) {
  
  library(ggplot2)
  library(reshape2)
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  count_matrix <- df[, count_cols, drop = FALSE]
  
  # Transform if requested
  if (log_transform) {
    count_matrix <- log2(count_matrix + 1)
    ylab <- "log2(Count + 1)"
  } else {
    ylab <- "Count"
  }
  
  # Reshape for plotting
  count_long <- melt(as.matrix(count_matrix))
  colnames(count_long) <- c("Gene", "Sample", "Expression")
  
  # Create plot
  p <- ggplot(count_long, aes(x = Sample, y = Expression, fill = Sample)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylab(ylab) +
    ggtitle("Expression Distribution Across Samples")
  
  return(p)
}
