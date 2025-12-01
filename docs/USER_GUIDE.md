# RNA-seq Preprocessing Pipeline - User Guide

## Table of Contents
1. [Introduction](#introduction)
2. [Quick Start](#quick-start)
3. [Detailed Workflow](#detailed-workflow)
4. [Function Reference](#function-reference)
5. [Troubleshooting](#troubleshooting)

## Introduction

This pipeline provides a complete workflow for preprocessing RNA-seq count data from FeatureCounts output. It handles:

- **Data loading** from standard FeatureCounts format
- **Quality control** with multiple visualization options
- **Chromosome filtering** to remove problematic regions
- **Gene annotation** converting ENSEMBL IDs to symbols
- **Duplicate handling** with multiple strategies
- **Expression filtering** to remove noise
- **Gene type filtering** to focus on protein-coding genes

## Quick Start

### 1. Installation

```r
# Install required packages
source("install_packages.R")
```

### 2. Prepare Your Data

Place your FeatureCounts output file in the `data/` directory:
```
data/featureCount.txt
```

### 3. Run the Pipeline

```r
# Run complete preprocessing
source("examples/full_preprocessing.R")
```

### 4. Check Results

Outputs will be saved in:
- `results/tables/` - Filtered count matrices
- `results/figures/` - QC plots

## Detailed Workflow

### Step 1: Load Data

```r
source("R/load_data.R")

# Load FeatureCounts output
df <- load_featurecounts("data/featureCount.txt")

# Extract count matrix only
counts <- get_count_matrix(df)

# Get gene metadata
metadata <- get_metadata(df)
```

**What this does:**
- Reads FeatureCounts tab-delimited output
- Shows dimensions and preview
- Separates counts from metadata

### Step 2: Filter Chromosomes

```r
source("R/filter_chromosomes.R")

# Remove Y chromosome, MT, and unplaced contigs
df_filtered <- filter_chromosomes(df, keep_pattern = "^[0-9X]+")

# Check chromosome distribution
chr_dist <- check_chromosome_distribution(df_filtered)
```

**Why filter chromosomes?**
- Y chromosome genes show sex-specific expression
- MT genes have different expression dynamics
- Unplaced contigs (KI, GL) are unreliable

**Typical filtering:**
- **Keep:** Autosomes (1-22) + X chromosome
- **Remove:** Y, MT, KI, GL chromosomes

### Step 3: Gene Annotation

```r
source("R/annotate_genes.R")

# Convert ENSEMBL to gene symbols
df_annotated <- annotate_genes(df_filtered, 
                               organism = "org.Hs.eg.db",
                               remove_na = TRUE)

# Add gene type information
df_annotated <- annotate_gene_types(df_annotated)
```

**Annotation details:**
- Uses Bioconductor annotation databases
- Removes genes without symbol matches
- Identifies gene types (protein-coding, lncRNA, etc.)

**Common organisms:**
- Human: `org.Hs.eg.db`
- Mouse: `org.Mm.eg.db`
- Rat: `org.Rn.eg.db`

### Step 4: Handle Duplicates

Multiple ENSEMBL IDs may map to the same gene symbol. Choose a strategy:

```r
source("R/handle_duplicates.R")

# Option 1: Random selection (fastest)
df_unique <- remove_duplicates(df_annotated, method = "random")

# Option 2: Average expression (conservative)
df_unique <- remove_duplicates(df_annotated, method = "average")

# Option 3: Highest expressing isoform
df_unique <- remove_duplicates(df_annotated, method = "highest")

# Find which genes are duplicated
duplicates <- find_duplicates(df_annotated)
```

**Method comparison:**

| Method | Speed | Use Case |
|--------|-------|----------|
| `random` | Fastest | Quick analysis, low duplicate % |
| `average` | Medium | Conservative, preserves information |
| `highest` | Medium | Focus on dominant isoform |

### Step 5: Filter Low Expression

```r
source("R/filter_expression.R")

# Method 1: Filter by mean expression
df_filtered <- filter_low_expression(df_unique, min_mean = 10)

# Method 2: Filter by counts in samples
df_filtered <- filter_low_expression(df_unique, 
                                     min_count = 5, 
                                     min_samples = 3)

# Get expression statistics
expr_stats <- get_expression_stats(df_filtered)
```

**Filtering rationale:**
- Removes genes with insufficient signal
- Reduces multiple testing burden
- Improves statistical power

**Recommended thresholds:**
- Mean expression > 10 counts
- OR count > 5 in ≥ 3 samples

### Step 6: Filter by Gene Type

```r
# Keep only protein-coding genes
df_protein <- filter_by_gene_type(df_filtered, 
                                  gene_type = "protein-coding")

# Or keep multiple types
df_multi <- filter_by_gene_type(df_filtered, 
                                gene_type = c("protein-coding", "lncRNA"))
```

**Common gene types:**
- `protein-coding` - Standard protein-coding genes (~20,000 in human)
- `lncRNA` - Long non-coding RNAs
- `miRNA` - microRNAs
- `pseudogene` - Pseudogenes

### Step 7: Quality Control

```r
source("R/qc_plots.R")

# Gene type composition
plot_gene_type_composition(df_filtered, 
                          output_dir = "results/figures")

# Sample correlation
cor_matrix <- plot_sample_correlation(df_filtered)

# PCA
pca_result <- plot_pca(df_filtered, 
                       color_by = sample_groups)

# Expression distribution
p <- plot_expression_distribution(df_filtered)

# Complete QC report
qc_results <- create_qc_report(df_filtered)
```

**QC Checks:**
1. **Gene type composition** - Should be mostly protein-coding
2. **Sample correlation** - Replicates should cluster (r > 0.9)
3. **PCA** - Check for batch effects or outliers
4. **Expression distribution** - Should be similar across samples

### Step 8: Save Results

```r
source("R/utils.R")

# Save count matrix
save_count_matrix(df_filtered, 
                 "results/tables/filtered_counts.csv",
                 include_metadata = TRUE)

# Get sample information
sample_info <- get_sample_info(df_filtered)
write.csv(sample_info, "results/tables/sample_info.csv")

# Print summary
print_processing_summary(original_df, final_df)
```

## Function Reference

### Data Loading

| Function | Purpose |
|----------|---------|
| `load_featurecounts()` | Load FeatureCounts output |
| `get_count_matrix()` | Extract count columns only |
| `get_metadata()` | Extract gene metadata |

### Filtering

| Function | Purpose |
|----------|---------|
| `filter_chromosomes()` | Remove unwanted chromosomes |
| `filter_low_expression()` | Remove lowly expressed genes |
| `filter_by_gene_type()` | Keep specific gene types |

### Annotation

| Function | Purpose |
|----------|---------|
| `annotate_genes()` | ENSEMBL → Symbol conversion |
| `annotate_gene_types()` | Add gene type info |
| `set_symbol_rownames()` | Use symbols as row names |

### Duplicate Handling

| Function | Purpose |
|----------|---------|
| `remove_duplicates()` | Handle duplicate symbols |
| `find_duplicates()` | Identify duplicated genes |
| `compare_duplicate_methods()` | Compare strategies |

### QC & Visualization

| Function | Purpose |
|----------|---------|
| `plot_gene_type_composition()` | Gene type distribution |
| `plot_sample_correlation()` | Sample correlation heatmap |
| `plot_pca()` | PCA plot |
| `plot_expression_distribution()` | Expression distributions |
| `create_qc_report()` | Complete QC report |

### Utilities

| Function | Purpose |
|----------|---------|
| `rename_samples()` | Rename sample columns |
| `save_count_matrix()` | Save with metadata |
| `load_count_matrix()` | Load saved matrix |
| `get_sample_info()` | Extract sample statistics |
| `print_processing_summary()` | Show processing stats |

## Troubleshooting

### Problem: File not found

```
Error: File not found: data/featureCount.txt
```

**Solution:** Ensure your FeatureCounts output is in the correct location.

### Problem: Annotation fails

```
Error: Package org.Hs.eg.db is required
```

**Solution:** Install Bioconductor packages:
```r
BiocManager::install("org.Hs.eg.db")
```

### Problem: Too many duplicates

```
Duplicate entries: 5000
```

**Solution:** This is normal. Different ENSEMBL IDs can map to the same symbol. Use `remove_duplicates()` to handle them.

### Problem: Low gene counts after filtering

**Solution:** 
- Check filtering thresholds (may be too stringent)
- Verify input data quality
- Review QC plots for issues

### Problem: Memory issues with large datasets

**Solution:**
- Process samples in batches
- Use `method = "random"` for faster duplicate handling
- Increase R memory limit: `options(java.parameters = "-Xmx8g")`

## Best Practices

1. **Always save intermediate results** - Helps with debugging
2. **Generate QC plots** - Visual inspection catches issues
3. **Document your choices** - Record thresholds and methods
4. **Version control** - Use git to track changes
5. **Reproducibility** - Set random seeds: `set.seed(123)`

## Next Steps

After preprocessing:
1. **Differential Expression** - Use DESeq2 or edgeR
2. **Pathway Analysis** - GSEA, GO enrichment
3. **Visualization** - Volcano plots, heatmaps
4. **Integration** - Multi-omics analysis

## Citation

If you use this pipeline, please cite:

```
[Your Name] (2024). RNA-seq Preprocessing Pipeline. 
GitHub: https://github.com/yourusername/rnaseq-preprocessing
```

## Support

For issues or questions:
- GitHub Issues: https://github.com/yourusername/rnaseq-preprocessing/issues
- Email: your.email@example.com
