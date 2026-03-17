#' Compute Weighted Information Gain Using Bootstrapped Resampling
#'
#' Run IG_weighted_resample_fast_repro : computes weighted information gain (IG) for 
#' each gene using multiple bootstrap resamples of the dataset. For every 
#' bootstrap iteration, the function recomputes IG using class-specific weights 
#' (optional), and aggregates the results across all iterations using both the 
#' median and the mean.
#'
#' The function is optimized for speed using parallel processing and returns a 
#' gene-level summary of weighted IG estimates.
#' 
#' Helper function (not for users)
#'
#' @param df A data frame or matrix where rows are samples and columns are genes, 
#'   with one column representing the class label.
#' @param class_col Character. The name of the column in `df` that contains the 
#'   class labels.
#' @param B Integer. Number of bootstrap resamples to run (default: 300).
#' @param N_per_boot Integer. Number of samples to include in each bootstrap 
#'   resample. Defaults to `nrow(df)`.
#' @param class_weights Optional numeric vector giving weights for each class. 
#'   If `NULL`, classes are treated equally.
#' @param seed Integer for random number generation (default: 42) for reproducible 
#'   resampling.
#' @param n_cores_boot Integer. Number of CPU cores to use for parallel bootstrap 
#'   execution. Defaults to `detectCores() - 1` with a minimum of 1.
#'
#' @return A data frame with one row per gene containing: Gene.Symbol, Median weighted IG across bootstrap iterations,
#'  and Mean weighted IG across bootstrap iterations. 
#' @noRd
IG_weighted_resample_fast_repro <- function(df,
                                            class_col,
                                            B = 300,
                                            N_per_boot = nrow(df),
                                            class_weights = NULL,
                                            seed = 42,
                                            n_cores_boot = max(1L, parallel::detectCores(logical = TRUE) - 1)) {
  stopifnot(class_col %in% names(df))
  
  # reproducible RNG streams
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  boot_seeds <- sample.int(.Machine$integer.max, B)
  
  # class factor
  y <- as.factor(df[[class_col]])
  df[[class_col]] <- y
  feats <- setdiff(names(df), class_col)
  df_small <- df[, c(feats, class_col), drop = FALSE]
  N <- nrow(df_small)
  
  # weights = N/(K*N_y)
  if (is.null(class_weights)) {
    tabs <- table(y)
    class_weights <- (sum(tabs) / (length(tabs) * tabs))
    class_weights <- setNames(as.numeric(class_weights), names(tabs))
  }
  p <- class_weights[as.character(y)]
  p <- p / sum(p)
  
  form <- as.formula(paste(class_col, "~ ."))
  
  # worker
  .one_boot <- function(b) {
    set.seed(boot_seeds[b])
    idx <- sample.int(N, size = N_per_boot, replace = TRUE, prob = p)
    sub <- df_small[idx, , drop = FALSE]
    ig_df <- FSelectorRcpp::information_gain(form, data = sub, type = "infogain", threads = 1)
    ig_vec <- setNames(ig_df$importance, ig_df$attributes)
    ig_vec[feats]  # return in fixed order
  }
  
  # parallel but reproducible
  res_list <- parallel::mclapply(seq_len(B), .one_boot,
                                 mc.cores = n_cores_boot,
                                 mc.preschedule = FALSE,
                                 mc.set.seed = FALSE)
  
  ig_mat <- do.call(cbind, res_list)
  rownames(ig_mat) <- feats
  
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    IG_wrs_median <- matrixStats::rowMedians(ig_mat, na.rm = TRUE)
    IG_wrs_mean   <- matrixStats::rowMeans2(ig_mat, na.rm = TRUE)
  } else {
    IG_wrs_median <- apply(ig_mat, 1, median, na.rm = TRUE)
    IG_wrs_mean   <- rowMeans(ig_mat, na.rm = TRUE)
  }
  
  data.frame(
    Gene.Symbol     = feats,
    IG_wrs_median = IG_wrs_median,
    IG_wrs_mean   = IG_wrs_mean,
    row.names = NULL
  )
}
