#!/usr/bin/env Rscript
#
# Basic Tests for RNA-seq Preprocessing Functions
#

cat("Running basic tests...\n\n")

# Source all functions
source("R/load_data.R")
source("R/filter_chromosomes.R")
source("R/annotate_genes.R")
source("R/handle_duplicates.R")
source("R/filter_expression.R")
source("R/utils.R")

# Create test data
create_test_data <- function() {
  # Simple test count matrix
  test_counts <- data.frame(
    Chr = c("1", "2", "X", "Y", "MT"),
    Start = c(1000, 2000, 3000, 4000, 5000),
    End = c(2000, 3000, 4000, 5000, 6000),
    Strand = c("+", "-", "+", "-", "+"),
    Length = c(1000, 1000, 1000, 1000, 1000),
    Sample1 = c(100, 50, 200, 30, 500),
    Sample2 = c(110, 55, 210, 35, 510),
    Sample3 = c(105, 52, 205, 32, 505),
    row.names = c("ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG005")
  )
  return(test_counts)
}

# Test 1: Chromosome filtering
cat("Test 1: Chromosome filtering...\n")
test_df <- create_test_data()
filtered <- filter_chromosomes(test_df, keep_pattern = "^[0-9X]+")
stopifnot(nrow(filtered) == 3)  # Should keep 1, 2, X
cat("✓ Chromosome filtering works\n\n")

# Test 2: Count matrix extraction
cat("Test 2: Count matrix extraction...\n")
count_matrix <- get_count_matrix(test_df)
stopifnot(ncol(count_matrix) == 3)  # Should have 3 sample columns
stopifnot(all(c("Sample1", "Sample2", "Sample3") %in% colnames(count_matrix)))
cat("✓ Count matrix extraction works\n\n")

# Test 3: Expression filtering
cat("Test 3: Expression filtering...\n")
filtered_expr <- filter_low_expression(test_df, min_mean = 60)
stopifnot(nrow(filtered_expr) == 3)  # Should keep genes with mean > 60
cat("✓ Expression filtering works\n\n")

# Test 4: Sample info extraction
cat("Test 4: Sample information...\n")
sample_info <- get_sample_info(test_df)
stopifnot(nrow(sample_info) == 3)  # 3 samples
stopifnot("TotalCounts" %in% colnames(sample_info))
cat("✓ Sample info extraction works\n\n")

# Test 5: Duplicate handling with test data
cat("Test 5: Duplicate handling...\n")
dup_df <- data.frame(
  SYMBOL = c("GENE1", "GENE1", "GENE2"),
  Sample1 = c(100, 150, 200),
  Sample2 = c(110, 160, 210)
)
rownames(dup_df) <- c("ENSG001", "ENSG002", "ENSG003")

# Test random method
unique_random <- remove_duplicates(dup_df, method = "random", symbol_col = "SYMBOL")
stopifnot(nrow(unique_random) == 2)  # Should have 2 unique symbols
cat("✓ Duplicate handling works\n\n")

cat("=" , rep("=", 50), "\n", sep = "")
cat("All basic tests passed!\n")
cat("=" , rep("=", 50), "\n", sep = "")
