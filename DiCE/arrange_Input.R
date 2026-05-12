#' Convert expression data to DiCE format
#'
#' Converts expression and metadata tables into DiCE-compatible format,
#' where rows are samples/cells, columns are genes, and a `class`
#' column contains phenotype/treatment labels.
#'
#' @param metadata A metadata data frame containing sample/cell annotations.
#' @param exp_data A gene expression data frame or matrix.
#'
#' @return A DiCE-formatted data frame.
#'
#' @noRd
exp_to_DiCE_format <- function(metadata,exp_data){
  
  # Step 1: Detect sample/cell ID column in metadata
  sample_col_patterns <- c("sample_id", "sampleid", "sample", "cell_id", "cellid", "cell", "sample id")
  sample_col <- names(metadata)[tolower(names(metadata)) %in% sample_col_patterns][1]
  if (is.na(sample_col)) stop("No sample/cell ID column found in metadata.")
  
  # Step 2: Detect phenotype column in metadata
  phenotype_col_patterns <- c("phenotype", "treatment", "class", "group", "condition")
  phenotype_col <- names(metadata)[tolower(names(metadata)) %in% phenotype_col_patterns][1]
  if (is.na(phenotype_col)) stop("No phenotype/treatment column found in metadata.")
  
  # Step 3: Detect gene column in exp_data
  gene_col_patterns <- c("gene", "genes", "gene_name", "gene_names", "genename", "genenames",
                         "gene.name", "gene.names", "gene_symbol", "genesymbol", "gene.symbol",
                         "gene_symbols", "genesymbols", "gene.symbols", "symbol", "symbols")
  gene_col <- names(exp_data)[tolower(names(exp_data)) %in% gene_col_patterns][1]
  
  # Step 4: Transpose exp_data so rows = samples/cells, cols = genes
  if (!is.na(gene_col)) {
    rownames(exp_data) <- exp_data[[gene_col]]
    exp_data[[gene_col]] <- NULL
  }
  exp_data <- as.data.frame(t(exp_data))
  
  # Step 5: Match samples/cells and add class column from metadata
  metadata_indexed <- metadata[match(rownames(exp_data), metadata[[sample_col]]), ]
  exp_data$class <- metadata_indexed[[phenotype_col]]
  
  if (any(is.na(exp_data$class))) {
    warning("Some samples/cells in expression data could not be matched in metadata.")
  }
  
  return(exp_data)
}



