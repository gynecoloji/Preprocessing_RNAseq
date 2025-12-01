#!/usr/bin/env Rscript
#
# Install Required R Packages for RNA-seq Preprocessing Pipeline
#

cat("Installing required R packages...\n\n")

# Function to install packages if not already installed
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  } else {
    cat(pkg, "already installed\n")
  }
}

# CRAN packages
cat("=== Installing CRAN packages ===\n")
cran_packages <- c(
  "stringr",
  "ggplot2",
  "dplyr",
  "pheatmap",
  "ggrepel",
  "ggpubr",
  "RColorBrewer",
  "reshape2",
  "devtools"
)

for (pkg in cran_packages) {
  install_if_missing(pkg, bioc = FALSE)
}

# Bioconductor packages
cat("\n=== Installing Bioconductor packages ===\n")
bioc_packages <- c(
  "DESeq2",
  "edgeR",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, bioc = TRUE)
}

# Optional: Mouse annotation (if analyzing mouse data)
cat("\n=== Optional packages ===\n")
cat("Installing mouse annotation database (optional)...\n")
install_if_missing("org.Mm.eg.db", bioc = TRUE)

# Verify installations
cat("\n=== Verifying installations ===\n")
all_packages <- c(cran_packages, bioc_packages)
success <- TRUE

for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "FAILED\n")
    success <- FALSE
  }
}

if (success) {
  cat("\n✓ All required packages installed successfully!\n")
  cat("\nYou can now run the preprocessing pipeline:\n")
  cat("  Rscript examples/full_preprocessing.R\n")
} else {
  cat("\n✗ Some packages failed to install. Please install them manually.\n")
}
