# ------- Phase I: Construction of a candidate gene pool by Differential Gene Expression -------
# analysis with a loose cutoff
#' Run Phase 1: Construction of a candidate gene pool by Differential Gene Expression analysis with a loose cutoff.
#' Helper function (not for users)
#'
#' @param dge_data Dataframe of DGE results
#' @param loose_criteria Column name for p-value or adj.p-value
#' @param loose_cutoff Numeric cutoff
#' @param logFC_cutoff Log fold change threshold
#'
#' @return Filtered dataframe of Phase1 candidate genes
#' @noRd
run_phase1 <- function(dge_data,loose_criteria,loose_cutoff,logFC_cutoff){
  filtered_df <- dge_data[dge_data[[loose_criteria]] <= loose_cutoff &
                            abs(dge_data$logFC) >= logFC_cutoff, ]
  filtered_df$Phase <- "I"
  return(filtered_df)
}