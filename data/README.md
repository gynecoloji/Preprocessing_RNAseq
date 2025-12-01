# Data Directory

Place your FeatureCounts output file here.

## Expected Format

The pipeline expects a FeatureCounts output file with the following structure:

```
Geneid  Chr     Start   End     Strand  Length  Sample1.bam  Sample2.bam  ...
ENSG00000000003 X       99883667        99894988        -       4390    156             203
ENSG00000000005 X       99839799        99854882        +       1610    0               1
...
```

## File Requirements

- **Format**: Tab-delimited text file
- **Header**: Must include column names
- **Geneid**: ENSEMBL gene IDs as first column
- **Metadata**: Chr, Start, End, Strand, Length columns
- **Counts**: One column per sample (BAM file)

## Example Files

If you don't have data yet, you can download example data:

```bash
# Example: Download from GEO or your data repository
```

## Generating FeatureCounts Output

If you're starting from fastq files, generate counts using my [Snakemake_RNAseq](https://github.com/gynecoloji/SnakeMake_RNAseq)
