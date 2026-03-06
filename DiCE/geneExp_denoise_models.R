# Get low-dimensional representations of single-cell RNAseq data

#' Run ZINB-WaVE : Zero-Inflated Negative Binomial - Weighted Adaptation for Variance and Effects Generalized linear model (GLM)
#' Helper function (not for users)
#'
#' @param raw_geneExp Raw UMI count matrix
#'
#' @return Dataframe of denoised gene expression values
#' @noRd
zinbWave_model <- function(raw_geneExp,min_count=1,min_cells=1){
  raw_geneExp_tmp <- raw_geneExp[,-ncol(raw_geneExp)]
  class <- raw_geneExp[,ncol(raw_geneExp)]
  
  count_data <- t(raw_geneExp_tmp)
  count_data <- round(as.matrix(count_data))
  
  # remove genes that have zero counts in all cells 
  keep_genes <- rowSums(count_data >= min_count) >= min_cells
  counts_filtered <- count_data[keep_genes, ]
  
  # create singlecellexperiment object
  # exp_obj <- SingleCellExperiment(assays = list(counts = counts_filtered))
  
  # wrap into SE
  exp_obj <- SummarizedExperiment(assays = list(counts = counts_filtered))
  
  # run ZINB-WaVE
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))
  
  zinb <- zinbwave(exp_obj,
                   K=0,
                   epsilon = 1e12,
                   normalizedValues = TRUE,
                   verbose = FALSE,
                   BPPARAM = BiocParallel::MulticoreParam(workers = 4))
  
  # extract denoised gene expression
  denoised <- assay(zinb,"normalizedValues") # gene x cells
  denoised <- log2(denoised + 1)
  
  denoised_t <- t(denoised)
  
  # Convert to data frames (if not already)
  denoised_t_df <- as.data.frame(denoised_t)
  
  # Extract the class column with matching rownames
  class_col <- colnames(raw_geneExp)[ncol(raw_geneExp)]
  matched_class <- raw_geneExp[rownames(denoised_t_df), class_col]
  
  # Add class as a new column, ensuring alignment
  denoised_t_df$class <- matched_class
  
  return(denoised_t_df)
}

#' Denoise scRNA-seq expression using NewWave and log-normalization
#' @param raw_geneExp A data.frame or matrix with cells as rows and
#'   genes as columns. The last column must contain class labels.
#' @param min_count Minimum count threshold to retain a gene.
#'   Genes must have counts >= \code{min_count} in at least
#'   \code{min_cells} cells.
#' @param min_cells Minimum number of cells required for a gene
#'   to be retained.
#' @param workers Number of parallel workers to use in NewWave.
#' @param K Number of latent factors to estimate in NewWave.
#'   Must be >= 1. Default is 5.
#'
#' @return A data.frame with Log-normalized gene expression values (cells × genes) and class column
#' @noRd
newWave_model <- function(raw_geneExp,
                          min_count = 1,
                          min_cells = 1,
                          workers = 4,
                          K = 5) {
  
  # split expression + class
  raw_geneExp_tmp <- raw_geneExp[, -ncol(raw_geneExp), drop = FALSE]
  class_col_name  <- colnames(raw_geneExp)[ncol(raw_geneExp)]
  
  # counts: genes x cells
  count_data <- t(raw_geneExp_tmp)
  count_data <- round(as.matrix(count_data))
  
  # filter genes
  keep_genes <- rowSums(count_data >= min_count) >= min_cells
  counts_filtered <- count_data[keep_genes, , drop = FALSE]
  
  # build SCE
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_filtered)
  )
  
  
  # run NewWave
  # K: number of latent factors (can be 0 if you only want covariate-based normalization) :contentReference[oaicite:1]{index=1}
  sce_nw <- NewWave::newWave(
    sce,
    K = K,
    verbose = FALSE,
    children = workers
  )
  
  # Use log-normalized counts for gene-gene correlation networks
  sce_nw <- scuttle::logNormCounts(sce_nw)
  mat <- SummarizedExperiment::assay(sce_nw, "logcounts")  # genes x cells
  out_df <- as.data.frame(t(mat))                          # cells x genes
  out_df$class <- raw_geneExp[rownames(out_df), class_col_name]
  
  return(out_df)
}


