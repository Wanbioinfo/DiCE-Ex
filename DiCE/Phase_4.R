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

#' Run DiCE Phase 4: Calculate centrality differences and ranks
#'
#' Helper function (not for users)
#'
#' @param centralities_df Data frame containing treatment and control
#'   centrality values (e.g., \code{"<Prefix>_treatment"}, \code{"<Prefix>_control"}).
#' @param centrality_list Character vector of centrality metrics to evaluate.
#'
#' @return A data frame containing the original input columns plus
#'   per-centrality \code{"<Prefix>_diff"} and \code{"<Prefix>_rank"} columns.
#'
#' @noRd
run_phase4 <- function(phase3_interactions_df, centrality_list) {
  
  treatment_interactions <- phase3_interactions_df[, c("source", "target", "weight_treatment")]
  colnames(treatment_interactions)[colnames(treatment_interactions) == "weight_treatment"] <- "weight"
  
  control_interactions <- phase3_interactions_df[, c("source", "target", "weight_control")]
  colnames(control_interactions)[colnames(control_interactions) == "weight_control"] <- "weight"
  
  # treatment centralities
  centralities_treatment <- calculate_centralities(treatment_interactions,centrality_list)
  centralities_treatment <- as.data.frame(centralities_treatment)
  suffix <- "_treatment" 
  colnames(centralities_treatment)[-1] <- paste0(colnames(centralities_treatment)[-1], suffix)
  
  # control centralities
  centralities_control <- calculate_centralities(control_interactions,centrality_list)
  centralities_control <- as.data.frame(centralities_control)
  suffix <- "_control" 
  colnames(centralities_control)[-1] <- paste0(colnames(centralities_control)[-1], suffix)
  
  
  # Merge treatment and control centralities
  centralities_df <- merge(centralities_treatment, centralities_control, by = "Gene.Symbol")
  
  # calculate treatment-control differences and centrality ranks
  centralities_df <- get_centDiff(centralities_df, centrality_list)
  centralities_df <- as.data.frame(centralities_df)
  
  # rearrange columns
  base_cols <- c(
    "Gene.Symbol", 
    if ("STRING_id" %in% colnames(centralities_df)) "STRING_id",
    if ("BioGRID_id" %in% colnames(centralities_df)) "BioGRID_id",
    "logFC", "adj.P.Val", "P.Value",
    if ("IG" %in% colnames(centralities_df)) "IG"
  )
  
  cent_cols <- unlist(lapply(centrality_list, function(cn) {
    pfx <- centrality_prefix(cn)
    c(
      paste0(pfx, "_treatment"),
      paste0(pfx, "_control"),
      paste0(pfx, "_diff"),
      paste0(pfx, "_rank")
    )
  }), use.names = FALSE)
  
  final_cols <- c(base_cols, cent_cols)
  final_cols <- intersect(final_cols, colnames(centralities_df))
  
  centralities_df <- centralities_df[, final_cols, drop = FALSE]
  
  return(centralities_df)
}