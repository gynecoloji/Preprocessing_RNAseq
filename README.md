# RNA-seq Data Preprocessing Pipeline

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A robust and modular R package for preprocessing RNA-seq count data from FeatureCounts output. Designed for computational biologists who need clean, reproducible data preparation workflows.

## Overview

This preprocessing pipeline handles the critical initial steps of RNA-seq analysis:
- **Data loading** from FeatureCounts output
- **Gene annotation** with ENSEMBL to Symbol conversion
- **Chromosome filtering** (remove Y, MT, and unplaced contigs)
- **Duplicate handling** with multiple strategies
- **Low-expression filtering** 
- **Gene type filtering** (protein-coding genes)
- **Quality control visualizations**

## Key Features

✅ **Modular Design** - Each preprocessing step is a separate, testable function  
✅ **Flexible Filtering** - Multiple strategies for duplicate genes and low-expression filtering  
✅ **Comprehensive QC** - Gene type composition plots and expression distribution analysis  
✅ **Well-Documented** - Clear function documentation and example workflows  
✅ **Production-Ready** - Error handling, logging, and reproducible outputs  

## Installation

### Prerequisites
- R ≥ 4.0
- Bioconductor packages

### Quick Start

```r
# Install required packages
source("install_packages.R")

# Or install manually
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "org.Hs.eg.db"))
install.packages(c("stringr", "ggplot2", "dplyr", "pheatmap", 
                   "ggrepel", "ggpubr", "RColorBrewer"))
```

## Project Structure

```
rnaseq-preprocessing/
├── R/                          # Core R functions
│   ├── load_data.R            # Data loading utilities
│   ├── filter_chromosomes.R   # Chromosome filtering
│   ├── annotate_genes.R       # Gene annotation functions
│   ├── handle_duplicates.R    # Duplicate gene strategies
│   ├── filter_expression.R    # Low-expression filtering
│   ├── qc_plots.R             # Quality control visualizations
│   └── utils.R                # Helper functions
├── examples/                   # Example workflows
│   ├── full_preprocessing.R   # Complete pipeline example
│   └── custom_workflow.R      # Custom analysis example
├── data/                       # Input data directory
├── results/                    # Output directory
│   ├── figures/               # QC plots
│   └── tables/                # Filtered count matrices
├── tests/                      # Unit tests
└── docs/                       # Documentation
```

## Quick Usage

### Basic Workflow

```r
# Source all functions
source("R/load_data.R")
source("R/filter_chromosomes.R")
source("R/annotate_genes.R")
source("R/handle_duplicates.R")
source("R/filter_expression.R")
source("R/qc_plots.R")

# Load data
counts <- load_featurecounts("data/featureCount.txt")

# Filter chromosomes (keep only autosomes and X)
counts_filtered <- filter_chromosomes(counts, keep_pattern = "^[0-9X]+")

# Annotate with gene symbols
counts_annotated <- annotate_genes(counts_filtered, 
                                   organism = "org.Hs.eg.db",
                                   id_type = "ENSEMBL")

# Handle duplicates (choose your strategy)
counts_unique <- remove_duplicates(counts_annotated, method = "random")
# OR
counts_averaged <- remove_duplicates(counts_annotated, method = "average")

# Filter low-expression genes
counts_final <- filter_low_expression(counts_unique, min_mean = 10)

# Generate QC plots
plot_gene_type_composition(counts_final, output_dir = "results/figures")

# Save filtered count matrix
write.csv(counts_final, "results/tables/filtered_counts.csv")
```

### Complete Pipeline (One Command)

```r
# Run the complete preprocessing pipeline
source("examples/full_preprocessing.R")
```

## Detailed Usage

### 1. Load FeatureCounts Output

```r
counts <- load_featurecounts("data/featureCount.txt")
# Returns: count matrix with ENSEMBL IDs as rownames
```

### 2. Filter Chromosomes

```r
# Keep only standard chromosomes (1-22, X)
counts_filtered <- filter_chromosomes(counts, keep_pattern = "^[0-9X]+")

# Or explicitly exclude specific chromosomes
counts_filtered <- filter_chromosomes(counts, 
                                      exclude_pattern = "^(Y|MT|KI|GL)")
```

### 3. Gene Annotation

```r
# Convert ENSEMBL to Symbol for human
counts_annotated <- annotate_genes(counts_filtered, 
                                   organism = "org.Hs.eg.db",
                                   id_type = "ENSEMBL",
                                   key_type = "SYMBOL")

# For mouse
counts_annotated <- annotate_genes(counts_filtered, 
                                   organism = "org.Mm.eg.db",
                                   id_type = "ENSEMBL")
```

### 4. Handle Duplicate Gene Symbols

Two strategies available:

```r
# Strategy 1: Random selection (faster)
counts_random <- remove_duplicates(counts_annotated, method = "random")

# Strategy 2: Average expression (more conservative)
counts_averaged <- remove_duplicates(counts_annotated, method = "average")
```

### 5. Filter Low-Expression Genes

```r
# Filter by mean expression across all samples
counts_filtered <- filter_low_expression(counts, min_mean = 10)

# Filter by minimum count in minimum number of samples
counts_filtered <- filter_low_expression(counts, 
                                         min_count = 5, 
                                         min_samples = 3)
```

### 6. Filter by Gene Type (Optional)

```r
# Keep only protein-coding genes
counts_protein <- filter_by_gene_type(counts, 
                                      gene_type = "protein-coding",
                                      organism = "org.Hs.eg.db")
```

### 7. Quality Control Plots

```r
# Gene type composition
plot_gene_type_composition(counts, 
                          output_dir = "results/figures",
                          output_prefix = "QC")

# Expression distribution
plot_expression_distribution(counts, 
                            output_dir = "results/figures")

# Sample correlation heatmap
plot_sample_correlation(counts, 
                       output_dir = "results/figures")
```

## Function Reference

### Core Functions

| Function | Purpose |
|----------|---------|
| `load_featurecounts()` | Load count data from FeatureCounts output |
| `filter_chromosomes()` | Remove unwanted chromosomes (Y, MT, etc.) |
| `annotate_genes()` | Convert ENSEMBL IDs to gene symbols |
| `remove_duplicates()` | Handle duplicate gene symbols |
| `filter_low_expression()` | Remove lowly expressed genes |
| `filter_by_gene_type()` | Keep specific gene types (e.g., protein-coding) |
| `plot_gene_type_composition()` | Visualize gene type distribution |
| `rename_samples()` | Rename sample columns to readable names |

## Example Dataset

The pipeline is demonstrated using PEO4 cell line RNA-seq data:

**Experimental Design:**
- BSA control (n=3)
- JAG1 treatment (n=3)  
- PEO4 with different treatments:
  - CBP (n=3)
  - IXZ (n=3)
  - Combination (n=3)
  - Control (n=3)

## Best Practices

1. **Always save intermediate outputs** - Helps with debugging and reproducibility
2. **Use version control** - Track changes to your analysis
3. **Document your choices** - Record filtering thresholds and duplicate handling methods
4. **Generate QC plots** - Visual inspection catches issues early
5. **Keep raw data separate** - Never overwrite original count files

## Output Files

The pipeline generates:
- **Filtered count matrices** (CSV format)
- **QC plots** (PDF format)
- **Gene annotation tables**
- **Processing logs**

## Testing

```r
# Run unit tests
source("tests/test_all.R")
```

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functions
4. Submit a pull request

## Citation

If you use this pipeline, please cite:

```
[Your Name] (2024). RNA-seq Data Preprocessing Pipeline. 
GitHub: https://github.com/yourusername/rnaseq-preprocessing
```

## License

MIT License - see [LICENSE](LICENSE) file

## Contact

- **Author**: Your Name
- **Email**: your.email@example.com
- **GitHub**: [@yourusername](https://github.com/yourusername)

## References

- Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology*
- Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." *Bioinformatics*
