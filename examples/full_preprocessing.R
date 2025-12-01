#!/usr/bin/env Rscript
#
# Complete RNA-seq Preprocessing Pipeline
# 
# This script demonstrates the full preprocessing workflow from
# FeatureCounts output to filtered, annotated count matrices
#
# Author: Your Name
# Date: 2024

# Load required libraries ======================================================
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(dplyr)

# Source all functions =========================================================
script_dir <- "R"
source(file.path(script_dir, "load_data.R"))
source(file.path(script_dir, "filter_chromosomes.R"))
source(file.path(script_dir, "annotate_genes.R"))
source(file.path(script_dir, "handle_duplicates.R"))
source(file.path(script_dir, "filter_expression.R"))
source(file.path(script_dir, "qc_plots.R"))
source(file.path(script_dir, "utils.R"))

# Set paths ====================================================================
data_dir <- "data"
results_dir <- "results"
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")

# Create directories
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# STEP 1: Load FeatureCounts data ============================================
cat("\n=== STEP 1: Loading Data ===\n")

# Adjust this path to your FeatureCounts output file
input_file <- file.path(data_dir, "featureCount.txt")

# Check if file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file, 
             "\nPlease place your featureCount.txt in the data/ directory"))
}

# Load the data
df_raw <- load_featurecounts(input_file)

# Keep original for comparison later
df_original <- df_raw

# STEP 2: Filter chromosomes =================================================
cat("\n=== STEP 2: Filtering Chromosomes ===\n")

# Remove genes on Y chromosome and unplaced contigs (KI, GL, MT)
# Keep only autosomes (1-22) and X chromosome
df <- filter_chromosomes(df_raw, keep_pattern = "^[0-9X]+")

# STEP 3: Annotate genes with symbols ========================================
cat("\n=== STEP 3: Gene Annotation ===\n")

# Convert ENSEMBL IDs to gene symbols
df <- annotate_genes(df, 
                     organism = "org.Hs.eg.db",
                     id_type = "ENSEMBL",
                     key_type = "SYMBOL",
                     remove_na = TRUE)

# Check for duplicates
cat("\nChecking for duplicate gene symbols...\n")
duplicate_summary <- find_duplicates(df, symbol_col = "SYMBOL")

# STEP 4: Handle duplicate gene symbols ======================================
cat("\n=== STEP 4: Handling Duplicates ===\n")

# Strategy 1: Random selection (faster, used in original code)
cat("\nStrategy 1: Random selection\n")
df_random <- remove_duplicates(df, method = "random", symbol_col = "SYMBOL")

# Strategy 2: Average expression (more conservative)
cat("\nStrategy 2: Average expression\n")
df_average <- remove_duplicates(df, method = "average", symbol_col = "SYMBOL")

# For this pipeline, we'll continue with both strategies

# STEP 5: Filter low-expression genes ========================================
cat("\n=== STEP 5: Filtering Low Expression ===\n")

# Filter genes with mean expression < 10
df_random_filtered <- filter_low_expression(df_random, min_mean = 10)
df_average_filtered <- filter_low_expression(df_average, min_mean = 10)

# STEP 6: Annotate and filter by gene type (optional) =======================
cat("\n=== STEP 6: Gene Type Annotation ===\n")

# Add gene type information
df_random_filtered <- annotate_gene_types(df_random_filtered, 
                                          organism = "org.Hs.eg.db")
df_average_filtered <- annotate_gene_types(df_average_filtered, 
                                           organism = "org.Hs.eg.db")

# Generate gene type composition plots
cat("\nGenerating gene type composition plots...\n")
plot_gene_type_composition(df_random_filtered, 
                          output_dir = figures_dir,
                          output_prefix = "random_method")

plot_gene_type_composition(df_average_filtered, 
                          output_dir = figures_dir,
                          output_prefix = "average_method")

# Optional: Keep only protein-coding genes
cat("\nFiltering for protein-coding genes...\n")
df_random_protein <- filter_by_gene_type(df_random_filtered, 
                                         gene_type = "protein-coding",
                                         add_genetype = FALSE)

df_average_protein <- filter_by_gene_type(df_average_filtered, 
                                          gene_type = "protein-coding",
                                          add_genetype = FALSE)

# STEP 7: Rename samples to readable names ===================================
cat("\n=== STEP 7: Renaming Samples ===\n")

# Define sample name mapping (adjust based on your actual column names)
# Original BAM file names from your analysis
old_names <- c(
  # You'll need to check the actual column names in your file
  # These are placeholders - adjust to match your actual column names
)

new_names <- c(
  'BSA_R1', 'BSA_R2', 'BSA_R3',
  'JAG1_R1', 'JAG1_R2', 'JAG1_R3',
  'PEO4_CBP_R1', 'PEO4_CBP_R2', 'PEO4_CBP_R3',
  'PEO4_Combi_R1', 'PEO4_Combi_R2', 'PEO4_Combi_R3',
  'PEO4_Ctrl_R1', 'PEO4_Ctrl_R2', 'PEO4_Ctrl_R3',
  'PEO4_IXZ_R1', 'PEO4_IXZ_R2', 'PEO4_IXZ_R3'
)

# Get current count column names
metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "SYMBOL", "GENETYPE")
current_count_cols <- setdiff(colnames(df_random_protein), metadata_cols)

cat("\nCurrent sample names:\n")
cat(paste(current_count_cols, collapse = "\n"), "\n")

cat("\nTo rename samples, update the 'old_names' vector in this script\n")
cat("with your actual BAM file column names, then uncomment the renaming code.\n\n")

# Uncomment and modify these lines once you've set up old_names correctly:
# df_random_protein <- rename_samples(df_random_protein, 
#                                     old_names = old_names, 
#                                     new_names = new_names)
# df_average_protein <- rename_samples(df_average_protein, 
#                                      old_names = old_names, 
#                                      new_names = new_names)

# STEP 8: Remove metadata columns for final count matrix ====================
cat("\n=== STEP 8: Preparing Final Count Matrices ===\n")

# Remove metadata columns, keep only counts
final_random <- df_random_protein[, setdiff(colnames(df_random_protein), 
                                            c("Chr", "Start", "End", "Strand", 
                                              "Length", "GENETYPE"))]

final_average <- df_average_protein[, setdiff(colnames(df_average_protein), 
                                              c("Chr", "Start", "End", "Strand", 
                                                "Length", "GENETYPE"))]

# Keep SYMBOL as column for reference, but genes are already row names
final_random$SYMBOL <- rownames(final_random)
final_average$SYMBOL <- rownames(final_average)

# Reorder columns (SYMBOL first, then counts)
count_cols <- setdiff(colnames(final_random), "SYMBOL")
final_random <- final_random[, c("SYMBOL", count_cols)]
final_average <- final_average[, c("SYMBOL", count_cols)]

# STEP 9: Save filtered count matrices =======================================
cat("\n=== STEP 9: Saving Results ===\n")

# Save both versions (corresponding to df3 and df4 from original code)
save_count_matrix(final_random, 
                 file.path(tables_dir, "filtered_counts_random.csv"),
                 include_metadata = TRUE)

save_count_matrix(final_average, 
                 file.path(tables_dir, "filtered_counts_average.csv"),
                 include_metadata = TRUE)

# STEP 10: Generate QC report ===============================================
cat("\n=== STEP 10: Quality Control ===\n")

# Create comprehensive QC plots
qc_results <- create_qc_report(df_random_filtered, 
                               output_dir = figures_dir)

# Print processing summary
print_processing_summary(df_original, final_random)

# Print sample information
cat("\n=== Sample Information ===\n")
sample_info <- get_sample_info(final_random)
print(sample_info)

# Save sample information
write.csv(sample_info, 
         file.path(tables_dir, "sample_info.csv"),
         row.names = FALSE)

# DONE =======================================================================
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("PREPROCESSING COMPLETE!\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("\nOutput files:\n")
cat("  - Filtered count matrices:", tables_dir, "\n")
cat("  - QC plots:", figures_dir, "\n")
cat("  - Sample information:", file.path(tables_dir, "sample_info.csv"), "\n")
cat("\nNext steps:\n")
cat("  1. Review QC plots in", figures_dir, "\n")
cat("  2. Choose between random or average duplicate handling method\n")
cat("  3. Proceed to differential expression analysis\n")
cat("\n")
