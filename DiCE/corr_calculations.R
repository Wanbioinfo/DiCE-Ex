# Set of correlation calculation methods

#' Run filtered_cor: Custom correlation function excluding shared-zero cells
#' Helper function (not for users)
#'
#' @param x Vector of gene expression of gene1 in all samples
#' @param y Vector of gene expression of gene2 in all samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#' @param corr_pval_cutoff P-value cutoff for filtering correlations
#'
#' @return Correlation between two genes
#' @noRd
filtered_cor <- function(x, y, corr_method, corr_pval_cutoff) {
  
  keep <- !(x == 0 & y == 0)
  
  x_keep <- x[keep]
  y_keep <- y[keep]
  
  if (length(x_keep) < 3) return(NA_real_)
  
  # zero variance check
  if (sd(x_keep) == 0 || sd(y_keep) == 0) {
    return(0)
  }
  
  mat <- cbind(x_keep, y_keep)
  
  res <- Hmisc::rcorr(mat, type = corr_method)
  
  corr_val <- res$r[1,2]
  p_val <- res$P[1,2]
  
  if (is.na(corr_val)) corr_val <- 0
  if (is.na(p_val)) p_val <- 1
  
  # set to 0 if p-value > cutoff
  if (p_val > corr_pval_cutoff) corr_val <- 0
  
  return(corr_val)
}



#' Run traditional_corr: Pearson OR spearman correlation taking all cells
#' Helper function (not for users)
#'
#' @param class_geneExp Dataframe of gene expression of all genes in a specific group of samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#' @param corr_pval_cutoff P-value cutoff for filtering correlations
#'
#' @return Correlation matrix between genes
#' @noRd
traditional_corr <- function(class_geneExp, corr_method, corr_pval_cutoff) {
  
  res <- Hmisc::rcorr(as.matrix(class_geneExp), type = corr_method)
  
  cor_matrix <- res$r
  p_matrix <- res$P
  
  if (any(is.na(cor_matrix))) {
    warning("Some correlations were undefined due to zero variance; these have been set to 0.",
            call. = FALSE)
  }
  
  cor_matrix[is.na(cor_matrix)] <- 0
  p_matrix[is.na(p_matrix)] <- 1
  
  # set to 0 if p-value > cutoff
  cor_matrix[p_matrix > corr_pval_cutoff] <- 0
  
  diag(cor_matrix) <- 1
  
  return(cor_matrix)
}

#' Run remove_Zerocells_corr: Pearson OR spearman correlation removing cells
#' where gene expression is zero for both genes in a gene-pair
#' Helper function (not for users)
#'
#' @param class_geneExp Dataframe of gene expression of all genes in a specific group of samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#' @param corr_pval_cutoff P-value cutoff for filtering correlations
#'
#' @return Absolute correlation matrix between genes
#' @noRd
remove_Zerocells_corr <- function(class_geneExp, corr_method, corr_pval_cutoff) {
  n_cores <- 4
  
  gene_names <- colnames(class_geneExp)
  n <- length(gene_names)
  
  cor_matrix <- matrix(NA_real_, n, n)
  colnames(cor_matrix) <- rownames(cor_matrix) <- gene_names
  
  gene_pairs <- combn(n, 2, simplify = FALSE)
  
  compute_pair <- function(pair) {
    i <- pair[1]
    j <- pair[2]
    
    val <- tryCatch(
      filtered_cor(
        x = class_geneExp[, i],
        y = class_geneExp[, j],
        corr_method = corr_method,
        corr_pval_cutoff = corr_pval_cutoff
      ),
      error = function(e) 0
    )
    
    list(i = i, j = j, value = val)
  }
  
  results <- parallel::mclapply(gene_pairs, compute_pair, mc.cores = n_cores)
  
  for (res in results) {
    cor_matrix[res$i, res$j] <- res$value
    cor_matrix[res$j, res$i] <- res$value
  }
  
  diag(cor_matrix) <- 1
  cor_matrix[is.na(cor_matrix)] <- 0
  
  return(cor_matrix)
}

