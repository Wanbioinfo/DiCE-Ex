# ------------------------------- Create PPI network for Phase 3 -----------------------------------

#' Run Get PPIs: Extract the PPI from StringDB for Phase2 selected genes
#' Helper function (not for users)
#'
#' @param string_protInfo_file Filepath of the StringDB protein information file
#' @param string_ppi_file Filepath of the StringDB PPI file
#' @param phase2_res Dataframe of Phase 2 results
#'
#' @return Dataframe of PPIs
#' @return Dataframe of vertices with the Gene.Symbol and String_ID
#' @noRd
create_PPI_fromPhase2 <- function(string_protInfo_file, string_ppi_file, phase2_res, score_threshold = 400){
  
  # Interested protein list from phase 2
  phase2_proteins <- phase2_res$Gene.Symbol
  
  # Read protein info
  # ppi_info_df <- read.delim(string_protInfo_file, header = TRUE, stringsAsFactors = FALSE)
  ppi_info_df <- read.delim(string_protInfo_file,
                            sep = "\t",
                            quote = "",
                            header = TRUE,
                            stringsAsFactors = FALSE)
  
  # Map STRING IDs to gene symbols
  colnames(ppi_info_df)[1] <- "STRING_id"
  colnames(ppi_info_df)[2] <- "Gene.Symbol"
  mapped_proteins <- ppi_info_df[, c("STRING_id", "Gene.Symbol")]
  colnames(mapped_proteins) <- c("STRING_id", "Gene.Symbol")
  
  # Take protein information of phase2 proteins
  mapped_proteins <- mapped_proteins %>% filter(Gene.Symbol %in% phase2_proteins)
  
  # Read STRING PPI file
  ppi_df <- read.table(string_ppi_file,
                       header = TRUE,
                       sep = "",          # any whitespace
                       stringsAsFactors = FALSE,
                       quote = "",
                       comment.char = "")
  # Ensure the type
  ppi_df$protein1 <- as.character(ppi_df$protein1)
  ppi_df$protein2 <- as.character(ppi_df$protein2)
  ppi_df$combined_score <- as.numeric(ppi_df$combined_score)
  
  # Filter the PPI with combined score >= 400
  ppi_df <- ppi_df %>% filter(combined_score >= score_threshold)
  
  # Map gene symbols -> STRING IDs
  string_ids <- mapped_proteins$STRING_id[mapped_proteins$Gene.Symbol %in% phase2_proteins]
  
  # Filter PPI file: keep only rows where both proteins are in your list
  interactions <- subset(ppi_df,
                         protein1 %in% string_ids & protein2 %in% string_ids)
  
  vertices <- unique(c(interactions$protein1,interactions$protein2))
  length(vertices)
  
  # Enforce order: make sure protein1 < protein2
  interactions$p1 <- pmin(interactions$protein1, interactions$protein2)
  interactions$p2 <- pmax(interactions$protein1, interactions$protein2)
  
  # Keep only relevant columns
  interactions <- interactions[, c("p1", "p2", "combined_score")]
  
  # Drop duplicates
  interactions <- interactions[!duplicated(interactions[, c("p1", "p2")]), ]
  
  # Map STRING IDs back to gene symbols for readability
  interactions <- interactions %>%
    merge(mapped_proteins, by.x = "p1", by.y = "STRING_id") %>%
    merge(mapped_proteins, by.x = "p2", by.y = "STRING_id") 
  
  colnames(interactions)[colnames(interactions) == "Gene.Symbol.x"] <- "Gene1"
  colnames(interactions)[colnames(interactions) == "Gene.Symbol.y"] <- "Gene2"
  
  interactions <- interactions[,c("Gene1","Gene2")]
  
  merged_df <- merge(phase2_res,
                     mapped_proteins,
                     by = "Gene.Symbol")

  final_df <- merged_df[, c("Gene.Symbol",  "STRING_id", "logFC", "adj.P.Val", "P.Value", if ("IG" %in% colnames(merged_df)) "IG")]
  
  return (list(interactions = interactions,
               mapped_proteins = final_df))
}

