# ----------------- Aggregate ranks in all phases to final output file -----------------

#' Normalize Centrality Metric Names
#'
#' Run norm_key() - converts an input string (typically a centrality metric name
#' provided by the user) into a normalized form suitable for matching. It
#' lowercases the string and removes all non-alphabetic characters to create a
#' standard key used for dictionary lookup.
#' 
#' Helper function (not for users)
#'
#' @param x Character string to normalize.
#'
#' @return A character string containing only lowercase alphabetic characters.
#'
#' @noRd
norm_key <- function(x) gsub("[^a-z]", "", tolower(x))


#' Mapping of User-Friendly Centrality Names to Data Frame Column Prefixes
#'
#' prefix_map provides a dictionary that maps flexible or commonly misspelled
#' user-input centrality metric names (e.g., "betweeness", "eig", "authority")
#' to the standardized column prefixes used internally in the DiCE output
#' data frame (e.g., "Betweenness", "EigenVector", "Authority").
#'
#' This enables users to write shorthand or approximate names while ensuring the
#' system consistently identifies the correct centrality metric columns.
#'
#' @format A named character vector where names are normalized keys and values are
#'   standardized column prefixes.
#'
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

#' Resolve Centrality Metric Names to Standardized Column Prefixes
#'
#' Run centrality_prefix() - takes a user-provided centrality metric name and returns
#' the corresponding standardized prefix used in the DiCE centrality result table.
#' The function normalizes the input name using norm_key(), looks up the
#' associated prefix in prefix_map, and throws an informative error if no mapping
#' exists.
#'
#' This allows the DiCE API to accept friendly, shorthand, or misspelled metric
#' names while still mapping correctly to the underlying data structure.
#' 
#' Helper function (not for users)
#'
#' @param x Character. A centrality metric name (e.g., `"eig"`, `"betweenness"`,
#'   `"pagerank"`, `"authority"`, `"harmonic"`).
#'
#' @return A character string giving the standardized column prefix (e.g.,
#'   `"EigenVector"`, `"Betweenness"`, `"PageRank"`).
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

#' Run Final ranking : Aggregate ranks in all phases to make a final rank (ProductOfRank)
#' Helper function (not for users)
#'
#' @param dge_data Dataframe of DGE results of all genes
#' @param phase1_res Dataframe of Phase1 results
#' @param phase2_res Dataframe of Phase2 results
#' @param phase4_centralities_df Dataframe of Phase4 centralities results
#' @param centrality_list List of the network centralities used in Phase4
#'
#' @return Dataframe of final ranking of each gene
#' @noRd
createFinalRanking <- function(dge_data, phase1_res, 
                               phase2_res, phase4_centralities_df,
                               centrality_list){

  # build the rank column names from the centrality list
  rank_cols <- paste0(sapply(centrality_list, centrality_prefix), "_rank")
  
  # add Phase3 genes and DiCE genes information for the final rank
  final_rank_df <- phase4_centralities_df
  
  # add Phase2 genes information for the final rank
  # add only the genes found in phase2
  if (!is.null(phase2_res) && nrow(phase2_res) > 0) {
    
    phase2_cols <- c(
      "Gene.Symbol",
      "logFC",
      "P.Value",
      "adj.P.Val",
      if ("IG" %in% colnames(phase2_res)) "IG",
      "Phase"
    )
    
    phase2_res <- phase2_res[, phase2_cols]
    
    only_phase2_rows <- phase2_res %>%
      as_tibble() %>%
      filter(!(Gene.Symbol %in% final_rank_df$Gene.Symbol)) %>%
      mutate(!!!setNames(rep(list(nrow(phase2_res)), length(rank_cols)), 
                         rank_cols))
    
    final_rank_df <- bind_rows(final_rank_df, only_phase2_rows)
    
  }
  
  
  # add Phase1 genes information for the final rank
  # add only the genes found in phase1
  phase1_cols <- c("Gene.Symbol","logFC","P.Value","adj.P.Val",(if ("IG" %in% colnames(phase1_res)) "IG"),"Phase" )
  phase1_res <- phase1_res[,phase1_cols]
  
  
  only_phase1_rows <- phase1_res %>%
    filter(!(Gene.Symbol %in% final_rank_df$Gene.Symbol)) %>%
    # dynamically add one rank column per centrality, all set to worst rank = nrow(phase1_res)
    mutate(!!!setNames(rep(list(nrow(phase1_res)), length(rank_cols)), rank_cols))

  
  final_rank_df <- bind_rows(final_rank_df, only_phase1_rows)
  
  # Add all expressed genes information for the final rank
  # Add only the remaining expressed genes
  dge_data <- dge_data[,c("Gene.Symbol","logFC","P.Value","adj.P.Val")]
  
  remaining_genes_rows <- dge_data %>%
    filter(!(Gene.Symbol %in% final_rank_df$Gene.Symbol)) %>%
    # dynamically add one rank column per centrality, all set to worst rank = nrow(dge_data)
    mutate(!!!setNames(rep(list(nrow(dge_data)), length(rank_cols)), rank_cols))

  
  final_rank_df <- bind_rows(final_rank_df, remaining_genes_rows)
  
  # find all columns ending with "_rank"
  rank_cols <- grep("_rank$", colnames(final_rank_df), value = TRUE)
  
  # Calculate final ensemble ranking 
  # compute consensus ranking - multiplies the ranks (lower product = better overall rank)
  # (M+1-z)/M
  final_rank_df$ProductofRank <- ensemble_rank(final_rank_df[, rank_cols, drop = FALSE],
                                                 method = "ProductOfRank")
  
  final_rank_df <- final_rank_df[order(final_rank_df[,"ProductofRank"],
                                       decreasing = TRUE ),]
  final_rank_df$Ensemble_Rank <- seq_len(nrow(final_rank_df))
  
  # Force Phase1, Phase2 and DGE ensemble ranks to the highest one in the group  
  final_rank_df$Ensemble_Rank[final_rank_df$Phase == "I"] <- nrow(phase1_res)
  
  if ("II" %in% final_rank_df$Phase) {
    final_rank_df$Ensemble_Rank[final_rank_df$Phase == "II"] <- nrow(phase2_res)
  }
  
  final_rank_df$Ensemble_Rank[is.na(final_rank_df$Phase)] <- nrow(dge_data)
  
  # Replace NA with "-" in Phase column
  final_rank_df$Phase[is.na(final_rank_df$Phase)] <- "-"
  
  
  return(final_rank_df)
}


