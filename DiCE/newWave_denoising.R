#' Denoise scRNA-seq expression using NewWave (zinbWave) and log-normalization
#'
#' This function filters the raw gene expression matrix based on differential gene expression (DGE) 
#' results and denoises the filtered matrix using the ZINB-WaVE model.
#'
#' This returns a log-normalized gene expression matrix suitable for downstream
#' gene-gene correlation network construction.
#'
#' 
#' @param dge_file Path to the '.rds' file containing the differential gene expression data. The columns must include: "Gene.Symbol", "logFC", "P.Value", "adj.P.Val".
#' @param rawGeneExp_file Path to the '.rds' file with raw gene expression data (cells x genes + label column).
#' @param loose_criteria Statistical significance metric used for initial gene filtering ("P.Value" or "adj.P.Val"). Default "adj.P.Val".
#' @param loose_cutoff Numeric threshold for filtering based on 'loose_criteria'. Genes with values ≤ this cutoff are retained. Default 0.05.
#' @param logFC_cutoff Minimum absolute log2 fold change threshold for retaining genes. Default 0.
#' @param min_count Minimum count threshold to retain a gene.
#'   Genes must have counts >= \code{min_count} in at least
#'   \code{min_cells} cells. Default is 1.
#' @param min_cells Minimum number of cells required for a gene
#'   to be retained. Default is 1.
#' @param workers Number of parallel workers to use in NewWave. Default is 4.
#' @param K Number of latent factors to estimate in NewWave.
#'   Must be >= 1. Default is 5.
#'
#' @return A data.frame with:
#' \itemize{
#'   \item Log-normalized gene expression values (cells × genes)
#'   \item A \code{class} column appended as the last column
#' }
#'
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom NewWave newWave
#' @importFrom scuttle logNormCounts
#' @importFrom SummarizedExperiment assay
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' denoised_df <- zinbWaVE_newModel_denoising(dge_file = "path/to/dge_data.rds",
#'                                           rawGeneExp_file = "path/to/raw_gene_exp.rds",
#'                                           loose_criteria = "adj.P.Val",
#'                                           loose_cutoff = 0.05,
#'                                           logFC_cutoff = 1)
#' }
#' 
#'
#' @export
zinbWaVE_newModel_denoising <- function(dge_file = NULL,
                          rawGeneExp_file = NULL,
                          loose_criteria = "adj.P.Val",
                          loose_cutoff = 0.05,
                          logFC_cutoff = 0,
                          min_count = 1,
                          min_cells = 1,
                          workers = 4,
                          K = 5) {
  
  
  if (dge_file == NULL){
    stop("Invalid or missing 'dge_file'")
  }
  
  if (rawGeneExp_file == NULL){
    stop("Invalid or missing 'rawGeneExp_file'")
  }
  

  # Load DGE data and gene expression data
  dge_data <- read_any(dge_file)
  rawGeneExp_data <- read_any(rawGeneExp_file)
  
  # Change column names to DiCE col names
  dge_data <- normalize_dge_cols(dge_data)
  
  
  # Filter genes and expression
  filtered_dge <- dge_data[dge_data[[loose_criteria]] <= loose_cutoff &
                             abs(dge_data$logFC) >= logFC_cutoff, ]
  
  filtered_genes <- filtered_dge$Gene.Symbol
  gene_cols <- colnames(rawGeneExp_data)[-ncol(rawGeneExp_data)]
  
  keep_cols <- intersect(gene_cols, filtered_genes)
  label_col <- colnames(rawGeneExp_data)[ncol(rawGeneExp_data)]
  
  rawGeneExp_filtered <- rawGeneExp_data[, c(keep_cols, label_col)]
  
  # ZinbWaVE denosing
  zinbWave_denoised <- newWave_model(rawGeneExp_filtered,
                                     min_count,
                                     min_cells,
                                     workers,
                                     K)

  
  return(zinbWave_denoised)
}






