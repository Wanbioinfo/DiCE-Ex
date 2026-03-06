# --------------------------------- Phase 4: ------------------------------------------------

#' Normalize a centrality name into a matching lookup key
#'
#' Converts a user-provided centrality metric name to a lowercase alphabetic key
#' by removing non-letters. Used internally for matching flexible input names.
#'
#' Helper function (not for users)
#' @param x Character string representing a centrality metric name.
#' @return A lowercase normalized key.
#' @noRd
norm_key <- function(x) gsub("[^a-z]", "", tolower(x))


#' Mapping of user-friendly centrality names to standardized column prefixes
#'
#' Named character vector linking normalized centrality keys to the canonical
#' centrality prefix names used in DiCE result tables.
#'  
#' For internal use 
#' @format A named character vector.
#' @noRd
prefix_map <- c(
  "betweeness"  = "Betweenness",   # common misspelling
  "betweenness" = "Betweenness",
  "eig"         = "EigenVector",
  "eigen vector" = "EigenVector",
  "eigenvector" = "EigenVector",
  "authority"   = "Authority",
  "strength"    = "Strength",
  "closeness"   = "Closeness",
  "pagerank"    = "PageRank",
  "harmonic"    = "Harmonic"
)

#' Resolve a user-provided centrality name to its standardized column prefix
#'
#' Matches a flexible or misspelled centrality name to the corresponding
#' standardized prefix using `norm_key()` and `prefix_map`.
#' 
#' Helper function (not for users)
#' 
#' @param x Character string specifying the centrality metric.
#' @return A standardized centrality column prefix.
#' 
#' @noRd
centrality_prefix <- function(x) {
  k <- norm_key(x)
  if (!k %in% names(prefix_map)) {
    stop("No prefix mapping for centrality: ", x,
         "\nAdd a mapping in prefix_map if your column prefix differs.")
  }
  unname(prefix_map[[k]])
}


#' Compute Treatment–Control Centrality Differences
#'
#' Calculates the difference between treatment and control centrality values for
#' each specified metric and assigns a rank based on the absolute difference.
#'
#' Helper function (not for users)
#' 
#' @param centralities_df A data frame containing centrality values with
#'   columns following the naming convention \code{"<Prefix>_treatment"} and
#'   \code{"<Prefix>_control"}.
#'
#' @param centrality_list Character vector of centrality metric names.
#'
#' @return The input data frame with additional \code{"<Prefix>_diff"} and
#'   \code{"<Prefix>_rank"} columns for each centrality metric.
#'
#' @noRd
get_centDiff <- function(centralities_df, centrality_list){
  
  for (centrality in centrality_list){
    centrality = tolower(centrality)
    
    pfx <- centrality_prefix(centrality)
    
    treat_col <- paste0(pfx,"_treatment")
    ctrl_col <- paste0(pfx,"_control")
    diff_col <- paste0(pfx, "_diff")
    rank_col <- paste0(pfx, "_rank")
    
    centralities_df[[diff_col]] <- centralities_df[[treat_col]] - centralities_df[[ctrl_col]]
    centralities_df <- centralities_df[order(abs(centralities_df[[diff_col]]), decreasing = TRUE), ]
    centralities_df[[rank_col]] <- seq_len(nrow(centralities_df))

  }
  
  return(centralities_df)
}


#' Compute a Numeric Cutoff Value
#'
#' Calculates a cutoff threshold from a numeric vector using one of three modes:
#' mean, median, or a top-K/top-K% percentile rule.
#'
#' Helper function (not for users)
#' 
#' @param x Numeric vector from which the cutoff is computed.
#' @param cutoff Method for selecting final candidate genes based on centrality metrics ("mean", "median", or "topK\%" such as "top10\%", "top25\%").
#'
#' @return A numeric cutoff threshold.
#'
#' @noRd
compute_cutoff <- function(x, cutoff) {
  x <- as.numeric(x)
  cutoff <- tolower(cutoff)
  if (cutoff == "mean") {
    return(mean(x, na.rm = TRUE))
  } else if (cutoff == "median") {
    return(median(x, na.rm = TRUE))
  } else if (grepl("^top[0-9]+%?$", cutoff)) {
    k <- as.numeric(gsub("top|%", "", cutoff))
    q <- 1 - (k / 100)  # keep top k% => threshold is (1 - k%) quantile
    return(quantile(x, probs = q, na.rm = TRUE, names = FALSE, type = 7))
  } else {
    stop("Invalid cutoff. Use 'mean', 'median', or 'topK'/'topK%'.")
  }
}


#' Identify DiCE-Passing Genes Based on Centrality Cutoffs
#'
#' Evaluates each gene across multiple centrality metrics by applying cutoff
#' thresholds to treatment and control centrality values, then selects genes
#' that pass a minimum number of metrics.
#'
#' Helper function (not for users)
#' 
#' @param centralities_df Data frame containing centrality values with
#'   \code{"<Prefix>_treatment"} and \code{"<Prefix>_control"} columns.
#' @param cutoff Cutoff rule passed to \code{compute_cutoff()} (e.g., "mean",
#'   "median", "top25\%").
#' @param centrality_list Character vector of centrality metrics to evaluate.
#' @param min_passCount Minimum number of centralities a gene must pass to be
#'   selected as a DiCE gene.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{pass_summary}: per-gene pass/fail and pass-count information
#'     \item \code{thresholds}: cutoff thresholds used for each centrality
#'     \item \code{DiCE_pass_df}: genes passing \code{min_passCount} metrics
#'   }
#'
#' @noRd
find_DiCEgenes <- function(centralities_df, cutoff, centrality_list, min_passCount){
  pass_flags <- list()
  thresholds  <- list() 
  
  # detect gene identifier column
  gene_col <- if ("Gene.Symbol" %in% names(centralities_df)) "Gene.Symbol" else
    if ("Gene" %in% names(centralities_df)) "Gene" else
      stop("Expected a column 'Gene.Symbol' or 'Gene'.")
  
  for (centrality in centrality_list){
    pfx <- centrality_prefix(centrality)
    treat_col <- paste0(pfx,"_treatment")
    ctrl_col <- paste0(pfx,"_control")
    
    # skip centralities not present
    if (!all(c(treat_col, ctrl_col) %in% colnames(centralities_df))) {
      warning(sprintf("Skipping %s (missing columns: %s or %s)", centrality, treat_col, ctrl_col))
      next
    }
    
    # thresholds
    treat_thr <- compute_cutoff(centralities_df[[treat_col]], cutoff)
    ctrl_thr <- compute_cutoff(centralities_df[[ctrl_col]], cutoff)
    thresholds[[pfx]] <- c(treatment = treat_thr, control = ctrl_thr)
    
    # pass if treatment >= treat_thr OR control >= ctrl_thr
    pass_flags[[pfx]] <- (centralities_df[[treat_col]] >= treat_thr) | (centralities_df[[ctrl_col]] >= ctrl_thr)
  }
  
  # combine results
  pass_df <- as.data.frame(pass_flags, optional = TRUE)
  pass_df[[gene_col]] <- centralities_df[[gene_col]]
  
  # pass count
  pass_df$pass_count <- rowSums(pass_df[names(pass_flags)], na.rm = TRUE)
  
  # list of passed centralities per gene
  pass_df$pass_centralities <- apply(
    pass_df[names(pass_flags)],
    1,
    function(x) {
      passed <- names(pass_flags)[which(x)]
      if (length(passed) == 0) return(NA_character_)
      paste(passed, collapse = ",")
    }
  )
  
  # DiCE genes = genes passing >= dice_min centralities
  DiCE_pass_df <- pass_df %>%
    filter(pass_count >= min_passCount)
  

  list(
    pass_summary = pass_df,        # includes per-gene pass info
    thresholds   = thresholds,     # per-centrality cutoffs used
    DiCE_pass_df   = DiCE_pass_df  # genes passing ≥ min_passCount
  )
  
}


#' Run DiCE Phase 4: Find DiCE genes based on network centrality measures
#'
#' Helper function (not for users)
#'
#' @param centralities_df Data frame containing treatment and control
#'   centrality values (e.g., \code{"<Prefix>_treatment"}, \code{"<Prefix>_control"}).
#' @param cutoff Cutoff rule used when identifying DiCE genes (passed to
#'   \code{compute_cutoff()}).
#' @param centrality_list Character vector of centrality metrics to evaluate.
#' @param min_passCount Minimum number of metrics a gene must pass to be
#'   classified as a DiCE gene.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{DiCE_genes_df}: Data frame of Phase 4–selected DiCE genes
#'     \item \code{phase4_centralities_df}: Full centrality table with ranks,
#'           ensemble ranking, and pass/fail annotations
#'   }
#'
#' @noRd
run_phase4 <- function(centralities_df, cutoff, centrality_list, min_passCount){
  
  # take difference between treatment and control centralities
  centralities_df <- get_centDiff(centralities_df, centrality_list)

  # find all columns ending with "_rank"
  rank_cols <- grep("_rank$", colnames(centralities_df), value = TRUE)
  
  # compute consensus ranking - multiplies the ranks (lower product = better overall rank)
  # (M+1-z)/M
  centralities_df$ProductofRank <- ensemble_rank(centralities_df[, rank_cols, drop = FALSE],
                                                 method = "ProductOfRank")
  
  # sort based on ProductofRank
  centralities_df <- centralities_df[order(centralities_df[,"ProductofRank"],
                                                         decreasing = TRUE ),]
  
  centralities_df$Phase4_rank <- seq_len(nrow(centralities_df))
  centralities_df <- as.data.frame(centralities_df)
  
  # find the DiCE genes
  dice_results <- find_DiCEgenes(centralities_df, cutoff, centrality_list, min_passCount)
  DiCE_pass_df <- dice_results$DiCE_pass_df
  pass_summary <- dice_results$pass_summary
  pass_summary <- pass_summary[,c("Gene.Symbol","pass_count","pass_centralities")]
  
  centralities_df <- merge(centralities_df, pass_summary, by = "Gene.Symbol")
  
  centralities_df <- centralities_df[order(centralities_df[,"ProductofRank"],
                                           decreasing = TRUE ),]
  
  # Rearrange the columns
  base_cols <- c("Gene.Symbol", "STRING_id", "logFC", "adj.P.Val", "P.Value", if ("IG" %in% colnames(centralities_df)) "IG")
  
  # for each requested centrality, ask for both treatment/control columns
  cent_cols <- unlist(lapply(centrality_list, function(cn) {
    pfx <- centrality_prefix(cn)
    c(paste0(pfx, "_treatment"), paste0(pfx, "_control"), 
      paste0(pfx, "_diff"),paste0(pfx,"_rank"))
  }), use.names = FALSE)
  
  tail_cols <- c("pass_count","pass_centralities","Phase")
  final_cols <- c(base_cols, cent_cols, tail_cols)
  
  centralities_df <- centralities_df[, final_cols, drop = FALSE]
  
  DiCE_genes <- DiCE_pass_df$Gene.Symbol
  centralities_df$Phase <- ifelse(centralities_df$Gene.Symbol %in% DiCE_genes, "DiCE", "III")
  
  DiCE_genes_df <- centralities_df[centralities_df$Gene.Symbol %in% DiCE_genes,]
  
  return(list(DiCE_genes_df = DiCE_genes_df,
              phase4_centralities_df = centralities_df))
  
}