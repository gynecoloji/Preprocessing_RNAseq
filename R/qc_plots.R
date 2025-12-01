#' Plot Gene Type Composition
#'
#' Visualize the proportion of different gene types across samples
#'
#' @param df Data frame with GENETYPE column and count data
#' @param output_dir Character. Directory to save plots (default: "results/figures")
#' @param output_prefix Character. Prefix for output files (default: "QC")
#' @param count_cols Character vector. Names of count columns (if NULL, auto-detect)
#' @return List of ggplot objects
#' @export

library(ggplot2)
library(stringr)
library(RColorBrewer)

plot_gene_type_composition <- function(df, 
                                       output_dir = "results/figures",
                                       output_prefix = "QC",
                                       count_cols = NULL) {
  
  if (!"GENETYPE" %in% colnames(df)) {
    stop("GENETYPE column not found. Run annotate_gene_types() first.")
  }
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  # Aggregate by gene type
  df_agg <- aggregate(. ~ GENETYPE, data = df[, c("GENETYPE", count_cols)], FUN = sum)
  
  # Convert to long format
  df_long <- reshape2::melt(df_agg, id.var = "GENETYPE")
  colnames(df_long) <- c("GENETYPE", "Sample", "Count")
  
  # Plot 1: Basic stacked bar
  p1 <- ggplot(df_long, aes(x = Sample, y = Count, fill = GENETYPE)) +
    geom_bar(stat = "identity", position = "fill") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Proportion") +
    ggtitle("Gene Type Composition Across Samples")
  
  ggsave(str_c(output_dir, "/", output_prefix, "_genetype_basic.pdf"), 
         p1, width = 8, height = 6)
  
  # Plot 2: Clean theme
  p2 <- ggplot(df_long, aes(x = Sample, y = Count, fill = GENETYPE)) +
    geom_bar(stat = "identity", position = "fill") +
    xlab("") +
    ylab("Percentage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank()) +
    ggtitle("Gene Type Distribution")
  
  ggsave(str_c(output_dir, "/", output_prefix, "_genetype_clean.pdf"), 
         p2, width = 8, height = 6)
  
  # Plot 3: With arrows
  p3 <- ggplot(df_long, aes(x = Sample, y = Count, fill = GENETYPE)) +
    geom_bar(stat = "identity", position = "fill") +
    xlab("") +
    ylab("Percentage") +
    scale_fill_brewer(palette = "Set3") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(color = "grey50", 
                                   arrow = arrow(type = "closed", 
                                                length = unit(0.1, 'inches')))) +
    ggtitle("Gene Type Composition")
  
  ggsave(str_c(output_dir, "/", output_prefix, "_genetype_styled.pdf"), 
         p3, width = 8, height = 6)
  
  cat("Saved 3 gene type composition plots to:", output_dir, "\n")
  
  return(list(basic = p1, clean = p2, styled = p3))
}


#' Plot Sample Correlation Heatmap
#'
#' Create correlation heatmap between samples
#'
#' @param df Data frame with count data
#' @param output_dir Character. Directory to save plot
#' @param output_file Character. Output filename (default: "sample_correlation.pdf")
#' @param count_cols Character vector. Names of count columns
#' @param method Character. Correlation method (default: "pearson")
#' @return Correlation matrix
#' @export

plot_sample_correlation <- function(df, 
                                   output_dir = "results/figures",
                                   output_file = "sample_correlation.pdf",
                                   count_cols = NULL,
                                   method = "pearson") {
  
  library(pheatmap)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  count_matrix <- df[, count_cols, drop = FALSE]
  
  # Log transform for correlation
  count_log <- log2(count_matrix + 1)
  
  # Calculate correlation
  cor_matrix <- cor(count_log, method = method)
  
  # Create heatmap
  pdf(file.path(output_dir, output_file), width = 10, height = 9)
  pheatmap(cor_matrix,
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = paste("Sample Correlation (", method, ")", sep = ""))
  dev.off()
  
  cat("Saved correlation heatmap to:", file.path(output_dir, output_file), "\n")
  
  return(cor_matrix)
}


#' Plot PCA
#'
#' Principal component analysis plot for sample clustering
#'
#' @param df Data frame with count data
#' @param output_dir Character. Directory to save plot
#' @param output_file Character. Output filename
#' @param count_cols Character vector. Names of count columns
#' @param color_by Character vector. Sample group labels for coloring
#' @return PCA results
#' @export

plot_pca <- function(df, 
                    output_dir = "results/figures",
                    output_file = "pca_plot.pdf",
                    count_cols = NULL,
                    color_by = NULL) {
  
  # Auto-detect count columns
  if (is.null(count_cols)) {
    metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
    count_cols <- setdiff(colnames(df), metadata_cols)
  }
  
  count_matrix <- df[, count_cols, drop = FALSE]
  
  # Log transform
  count_log <- log2(count_matrix + 1)
  
  # PCA
  pca_result <- prcomp(t(count_log), center = TRUE, scale. = TRUE)
  
  # Variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Create plot data
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = rownames(pca_result$x)
  )
  
  if (!is.null(color_by)) {
    pca_df$Group <- color_by
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
      geom_point(size = 4) +
      geom_text(vjust = -1, hjust = 0.5, size = 3)
  } else {
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(size = 4) +
      geom_text(vjust = -1, hjust = 0.5, size = 3)
  }
  
  p <- p +
    xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme_bw() +
    ggtitle("PCA of Samples")
  
  # Save
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ggsave(file.path(output_dir, output_file), p, width = 8, height = 6)
  
  cat("Saved PCA plot to:", file.path(output_dir, output_file), "\n")
  cat("PC1 variance explained:", round(var_explained[1], 1), "%\n")
  cat("PC2 variance explained:", round(var_explained[2], 1), "%\n")
  
  return(pca_result)
}


#' Create QC Report
#'
#' Generate comprehensive QC plots and summary statistics
#'
#' @param df Data frame with count data
#' @param output_dir Character. Directory to save outputs
#' @param sample_groups Character vector. Sample group labels (optional)
#' @return List with QC results
#' @export

create_qc_report <- function(df, 
                             output_dir = "results/figures",
                             sample_groups = NULL) {
  
  cat("Generating QC report...\n\n")
  
  results <- list()
  
  # Gene type composition (if GENETYPE available)
  if ("GENETYPE" %in% colnames(df)) {
    cat("1. Gene type composition plots\n")
    results$genetype_plots <- plot_gene_type_composition(df, output_dir = output_dir)
  }
  
  # Expression distribution
  cat("2. Expression distribution plot\n")
  p_expr <- plot_expression_distribution(df)
  ggsave(file.path(output_dir, "expression_distribution.pdf"), 
         p_expr, width = 10, height = 6)
  
  # Sample correlation
  cat("3. Sample correlation heatmap\n")
  results$correlation <- plot_sample_correlation(df, output_dir = output_dir)
  
  # PCA
  cat("4. PCA plot\n")
  results$pca <- plot_pca(df, output_dir = output_dir, color_by = sample_groups)
  
  cat("\nQC report complete! Check", output_dir, "for outputs.\n")
  
  return(results)
}
