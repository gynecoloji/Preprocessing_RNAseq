# Quick Reference Guide

## One-Line Commands

### Installation
```r
source("install_packages.R")
```

### Run Complete Pipeline
```r
source("examples/full_preprocessing.R")
```

### Run Tests
```r
source("tests/test_basic.R")
```

## Common Tasks

### Load Data
```r
source("R/load_data.R")
df <- load_featurecounts("data/featureCount.txt")
```

### Filter and Annotate
```r
# Filter chromosomes
source("R/filter_chromosomes.R")
df <- filter_chromosomes(df, keep_pattern = "^[0-9X]+")

# Annotate genes
source("R/annotate_genes.R")
df <- annotate_genes(df, organism = "org.Hs.eg.db")
```

### Handle Duplicates
```r
source("R/handle_duplicates.R")
df <- remove_duplicates(df, method = "random")  # or "average" or "highest"
```

### Filter Expression
```r
source("R/filter_expression.R")
df <- filter_low_expression(df, min_mean = 10)
```

### Generate QC Plots
```r
source("R/qc_plots.R")
create_qc_report(df, output_dir = "results/figures")
```

### Save Results
```r
source("R/utils.R")
save_count_matrix(df, "results/tables/filtered_counts.csv")
```

## File Locations

| What | Where |
|------|-------|
| Input data | `data/featureCount.txt` |
| Results | `results/tables/` |
| Figures | `results/figures/` |
| Functions | `R/` |
| Examples | `examples/` |

## Typical Workflow

```r
# 1. Load all functions
source("R/load_data.R")
source("R/filter_chromosomes.R")
source("R/annotate_genes.R")
source("R/handle_duplicates.R")
source("R/filter_expression.R")
source("R/qc_plots.R")
source("R/utils.R")

# 2. Process data
df <- load_featurecounts("data/featureCount.txt")
df <- filter_chromosomes(df, keep_pattern = "^[0-9X]+")
df <- annotate_genes(df, organism = "org.Hs.eg.db")
df <- remove_duplicates(df, method = "random")
df <- filter_low_expression(df, min_mean = 10)

# 3. QC and save
create_qc_report(df, output_dir = "results/figures")
save_count_matrix(df, "results/tables/filtered_counts.csv")
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Package not found | `source("install_packages.R")` |
| File not found | Check path: `data/featureCount.txt` |
| Memory error | Use `method = "random"` for duplicates |
| Too few genes | Lower filtering thresholds |

## Getting Help

- Read full guide: `docs/USER_GUIDE.md`
- Check examples: `examples/full_preprocessing.R`
- Open issue: GitHub Issues
- Email: your.email@example.com
