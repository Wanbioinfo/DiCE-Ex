#' Annotate transcription factors (internal use only)
#'
#' Internal helper function to annotate DiCE results with
#' high-confidence transcription factors (TFs) based on
#' UniProt annotation filters.
#'
#' Filters the provided TF reference file using:
#' \itemize{
#'   \item Reviewed entries
#'   \item Gene Ontology term: "DNA-binding transcription factor activity"
#'   \item DNA binding annotation
#'   \item Nuclear subcellular localization
#' }
#'
#' Adds a column \code{Is_TF} to the input results data frame.
#'
#' @param dice_results_df Data frame containing DiCE-ranked genes.
#' @param tf_file Path to UniProt-based TF annotation Excel file.
#'
#' @return Data frame with an additional \code{Is_TF} column.
#'
#' @noRd
annotateTFs <- function(dice_results_df, tf_file){
  
  df <- read_excel(tf_file)
  
  filtered_df <- df %>%
    filter(
      
      str_to_lower(Reviewed) == "reviewed",
      
      str_detect(str_to_lower(`Gene Ontology (GO)`),
                 "dna-binding transcription factor activity"),
      
      str_detect(str_to_lower(`DNA binding`), "dna_bind"),
      
      str_detect(str_to_lower(`Subcellular location [CC]`), "nucleus")
    )
  
  tf_genes <- filtered_df[["Gene Names (primary)"]]
  
  dice_results_df$Is_TF <- ifelse(
    dice_results_df$Gene.Symbol %in% tf_genes,
    "TF",
    NA
  )
  return(dice_results_df)
}