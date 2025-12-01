#' Annotate Genes with Symbols
#'
#' Convert ENSEMBL IDs to gene symbols using organism-specific annotation database
#'
#' @param df Data frame with ENSEMBL IDs as rownames
#' @param organism Character. Organism database (default: "org.Hs.eg.db" for human)
#' @param id_type Character. Input ID type (default: "ENSEMBL")
#' @param key_type Character. Output key type (default: "SYMBOL")
#' @param remove_na Logical. Remove genes without symbol annotation (default: TRUE)
#' @return Data frame with SYMBOL column added
#' @export
#' @examples
#' # Human genes
#' annotated <- annotate_genes(counts, organism = "org.Hs.eg.db")
#' 
#' # Mouse genes
#' annotated <- annotate_genes(counts, organism = "org.Mm.eg.db")

annotate_genes <- function(df, 
                          organism = "org.Hs.eg.db",
                          id_type = "ENSEMBL",
                          key_type = "SYMBOL",
                          remove_na = TRUE) {
  
  # Load organism database
  if (!requireNamespace(organism, quietly = TRUE)) {
    stop(paste("Package", organism, "is required. Install with: BiocManager::install('", 
               organism, "')", sep = ""))
  }
  
  library(organism, character.only = TRUE)
  
  cat("Annotating", nrow(df), "genes using", organism, "\n")
  
  # Get gene symbols
  annotation_genes <- AnnotationDbi::mapIds(
    get(organism),
    keys = rownames(df),
    column = key_type,
    keytype = id_type,
    multiVals = "first"  # Take first match if multiple
  )
  
  # Add to dataframe
  df$SYMBOL <- unname(annotation_genes)
  
  # Count matches
  n_matched <- sum(!is.na(df$SYMBOL))
  n_unmatched <- sum(is.na(df$SYMBOL))
  
  cat("Matched:", n_matched, "genes\n")
  cat("Unmatched:", n_unmatched, "genes\n")
  cat("Match rate:", round(n_matched / nrow(df) * 100, 2), "%\n")
  
  # Remove NA if requested
  if (remove_na) {
    n_before <- nrow(df)
    df <- df[complete.cases(df$SYMBOL), ]
    cat("Removed", n_before - nrow(df), "genes without symbols\n")
  }
  
  # Show unique symbols vs total rows
  n_unique <- length(unique(df$SYMBOL))
  cat("\nUnique gene symbols:", n_unique, "\n")
  cat("Total rows:", nrow(df), "\n")
  
  if (n_unique < nrow(df)) {
    n_duplicates <- nrow(df) - n_unique
    cat("Duplicate symbols:", n_duplicates, 
        "(will need to handle duplicates)\n")
  }
  
  return(df)
}


#' Annotate Gene Types
#'
#' Add gene type information (protein-coding, lncRNA, etc.)
#'
#' @param df Data frame with gene symbols as rownames or SYMBOL column
#' @param organism Character. Organism database (default: "org.Hs.eg.db")
#' @param symbol_col Character. Name of symbol column (default: "SYMBOL")
#' @return Data frame with GENETYPE column added
#' @export

annotate_gene_types <- function(df, 
                                organism = "org.Hs.eg.db",
                                symbol_col = "SYMBOL") {
  
  # Load organism database
  if (!requireNamespace(organism, quietly = TRUE)) {
    stop(paste("Package", organism, "is required"))
  }
  
  library(organism, character.only = TRUE)
  
  # Determine which symbols to annotate
  if (symbol_col %in% colnames(df)) {
    symbols <- df[[symbol_col]]
  } else {
    symbols <- rownames(df)
  }
  
  cat("Annotating gene types for", length(symbols), "genes\n")
  
  # Get gene types
  annotation_types <- AnnotationDbi::mapIds(
    get(organism),
    keys = symbols,
    column = "GENETYPE",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  df$GENETYPE <- unname(annotation_types)
  
  # Summary of gene types
  cat("\nGene type distribution:\n")
  print(table(df$GENETYPE, useNA = "ifany"))
  
  return(df)
}


#' Set Gene Symbols as Rownames
#'
#' Replace ENSEMBL IDs with gene symbols as rownames
#'
#' @param df Data frame with SYMBOL column
#' @param symbol_col Character. Name of symbol column (default: "SYMBOL")
#' @param keep_symbol_col Logical. Keep SYMBOL column (default: FALSE)
#' @return Data frame with symbols as rownames
#' @export

set_symbol_rownames <- function(df, 
                               symbol_col = "SYMBOL",
                               keep_symbol_col = FALSE) {
  
  if (!symbol_col %in% colnames(df)) {
    stop(paste("Column", symbol_col, "not found"))
  }
  
  # Check for duplicates
  if (any(duplicated(df[[symbol_col]]))) {
    warning("Duplicate symbols found. Consider using remove_duplicates() first.")
  }
  
  rownames(df) <- df[[symbol_col]]
  
  if (!keep_symbol_col) {
    df[[symbol_col]] <- NULL
  }
  
  return(df)
}
