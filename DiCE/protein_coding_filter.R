# Function to keep only protein coding genes

#' Retrieve protein-coding gene symbols for a given species
#'
#' Uses organism-specific annotation packages to return the set of
#' protein-coding genes based on the supplied NCBI taxonomy ID.
#' 
#' Helper function (not for users)
#'
#' @param taxonID Integer NCBI taxonomy ID (e.g., 9606 for human, 10090 for mouse).
#'
#' @return A character vector of unique gene symbols annotated as protein-coding.
#'
#' @noRd
get_protCodingGenes <- function(taxonID){
  
  # choose the correct annotation database object
  if (taxonID == 9606) {
    db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (taxonID == 10090) {
    db <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    stop("Unsupported taxonID: ", taxonID)
  }
  
  gene_info <- AnnotationDbi::select(
    db,
    keys    = AnnotationDbi::keys(db, keytype = "SYMBOL"),
    keytype = "SYMBOL",
    columns = c("GENENAME", "ENSEMBL", "GENETYPE")
  )
  
  protein_coding_genes <- gene_info %>%
    dplyr::filter(GENETYPE == "protein-coding") %>%
    dplyr::pull(SYMBOL) %>%
    unique()
  
  return(protein_coding_genes)
  
}


#' Subset DGE results to protein-coding genes
#'
#' Filters a DGE table to retain only rows whose \code{Gene.Symbol} corresponds
#' to protein-coding genes for the specified taxonomy ID.
#'
#' Helper function (not for users)
#' 
#' @param taxonID Integer NCBI taxonomy ID passed to \code{get_protCodingGenes()}.
#' @param dge_data Data frame of differential expression results with a
#'   \code{Gene.Symbol} column.
#'
#' @return A filtered DGE data frame containing only protein-coding genes.
#' @noRd
keep_protCoding_dgeData <- function(taxonID, dge_data){
  
  protein_coding_genes <- get_protCodingGenes(taxonID)
  dge_data <- dge_data %>% filter(Gene.Symbol %in% protein_coding_genes)
  
  return(dge_data)
  
}

#' Subset expression matrix to protein-coding genes
#'
#' Restricts an expression matrix to the specified protein-coding genes while
#' preserving the last column (e.g., class label or metadata).
#' 
#' Helper function (not for users)
#'
#' @param prot_genes Character vector of protein-coding gene symbols to retain.
#' @param exp_data Expression data frame or matrix with genes in columns and the
#'   last column typically representing a label or phenotype.
#'
#' @return An expression data frame containing only the selected genes and the
#'   final column.
#'
#' @noRd
keep_protCoding_expData <- function(prot_genes, exp_data){
  
  last_col <- names(exp_data)[ncol(exp_data)]
  final_cols <- c(prot_genes, last_col)
  exp_data <- exp_data[, names(exp_data) %in% final_cols, drop = FALSE]
  
  return(exp_data)
}

