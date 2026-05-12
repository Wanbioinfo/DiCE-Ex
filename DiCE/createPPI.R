# ------------------------------- Create PPI network for Phase 3 -----------------------------------

#' Run Get PPIs: Extract the PPI from StringDB for Phase2 selected genes
#' Helper function (not for users)
#'
#' @param string_protInfo_file Filepath of the StringDB protein information file
#' @param string_ppi_file Filepath of the StringDB PPI file
#' @param phase2_res Dataframe of Phase 2 results
#' @param stringDB_confidence Cutoff for confidence score in Stringdb PPI (Default is 400).
#' @return Dataframe of PPIs
#' @return Dataframe of vertices with the Gene.Symbol and String_ID
#' @noRd
create_stringPPI_fromPhase2 <- function(string_protInfo_file, 
                                  string_ppi_file, 
                                  stringDB_confidence = 400,
                                  phase2_res){
  
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
  ppi_df <- ppi_df %>% filter(combined_score >= stringDB_confidence)
  
  # Map gene symbols -> STRING IDs
  string_ids <- mapped_proteins$STRING_id[mapped_proteins$Gene.Symbol %in% phase2_proteins]
  
  # Filter PPI file: keep only rows where both proteins are in your list
  interactions <- subset(ppi_df,
                         protein1 %in% string_ids & protein2 %in% string_ids)
  
  vertices <- unique(c(interactions$protein1,interactions$protein2))
  
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
  
  final_df <- merged_df[, c("Gene.Symbol",  
                            if ("STRING_id" %in% colnames(merged_df)) "STRING_id",
                            if ("logFC" %in% colnames(merged_df)) "logFC", 
                            if ("adj.P.Val" %in% colnames(merged_df)) "adj.P.Val",
                            if ("P.Value" %in% colnames(merged_df)) "P.Value", 
                            if ("IG" %in% colnames(merged_df)) "IG")]
  
  return (list(interactions = interactions,
               mapped_proteins = final_df))
}

#' Run Get PPIs: Extract the PPI from BioGrid for Phase2 selected genes
#' Helper function (not for users)
#'
#' @param biogrid_ppi_file Filepath of the BioGrid PPI file
#' @param phase2_res Dataframe of Phase 2 results
#' @return Dataframe of PPIs
#' @return Dataframe of vertices with the Gene.Symbol and String_ID
#' @noRd
create_biogridPPI_fromPhase2 <- function(biogrid_ppi_file, phase2_res){
  
  # Interested protein list from phase 2
  phase2_proteins <- phase2_res$Gene.Symbol
  
  # Read BioGRID interaction file
  lines <- readLines(biogrid_ppi_file)
  header_line <- which(grepl("^INTERACTOR_A", lines))[1]
  
  biogrid_df <- read.delim(
    text = paste(lines[header_line:length(lines)], collapse = "\n"),
    sep = "\t",
    quote = "",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  
  # Ensure character types
  biogrid_df$OFFICIAL_SYMBOL_A <- as.character(biogrid_df$OFFICIAL_SYMBOL_A)
  biogrid_df$OFFICIAL_SYMBOL_B <- as.character(biogrid_df$OFFICIAL_SYMBOL_B)
  biogrid_df$ORGANISM_A_ID     <- as.integer(biogrid_df$ORGANISM_A_ID)
  biogrid_df$ORGANISM_B_ID     <- as.integer(biogrid_df$ORGANISM_B_ID)
  
  # Keep only human-human interactions (NCBI taxon ID 9606)
  biogrid_df <- biogrid_df %>%
    filter(ORGANISM_A_ID == 9606,
           ORGANISM_B_ID == 9606)
  
  
  # Filter: keep only rows where both interactors are in phase2 proteins
  interactions <- subset(biogrid_df,
                         OFFICIAL_SYMBOL_A %in% phase2_proteins &
                           OFFICIAL_SYMBOL_B %in% phase2_proteins)
  
  # Enforce canonical order to remove duplicates (A < B lexicographically)
  interactions$Gene1 <- pmin(interactions$OFFICIAL_SYMBOL_A,
                             interactions$OFFICIAL_SYMBOL_B)
  interactions$Gene2 <- pmax(interactions$OFFICIAL_SYMBOL_A,
                             interactions$OFFICIAL_SYMBOL_B)
  
  # Keep only gene pair columns
  interactions <- interactions[, c("Gene1", "Gene2")]
  
  # Remove self-loops
  interactions <- interactions %>%
    filter(Gene1 != Gene2)
  
  # Drop duplicate pairs
  interactions <- interactions[!duplicated(interactions[, c("Gene1", "Gene2")]), ]
  
  # Check number of unique proteins involved
  vertices <- unique(c(interactions$Gene1, interactions$Gene2))
  
  # Build a mapping of Gene Symbol -> BioGRID Interactor ID
  # Extract both directions (A and B) and combine
  biogrid_id_map <- rbind(
    biogrid_df[, c("INTERACTOR_A", "OFFICIAL_SYMBOL_A")] %>%
      rename(BioGRID_id = INTERACTOR_A, Gene.Symbol = OFFICIAL_SYMBOL_A),
    biogrid_df[, c("INTERACTOR_B", "OFFICIAL_SYMBOL_B")] %>%
      rename(BioGRID_id = INTERACTOR_B, Gene.Symbol = OFFICIAL_SYMBOL_B)
  ) %>%
    distinct(Gene.Symbol, .keep_all = TRUE)  # one BioGRID ID per gene symbol
  
  
  # Merge phase2 results with BioGRID IDs
  merged_df <- merge(phase2_res,
                     biogrid_id_map,
                     by = "Gene.Symbol")
  
  # Build final_df with BioGRID_id instead of STRING_id
  final_df <- merged_df[, c("Gene.Symbol",
                            if ("BioGRID_id"  %in% colnames(merged_df)) "BioGRID_id",
                            if ("logFC"       %in% colnames(merged_df)) "logFC",
                            if ("adj.P.Val"   %in% colnames(merged_df)) "adj.P.Val",
                            if ("P.Value"     %in% colnames(merged_df)) "P.Value",
                            if ("IG"          %in% colnames(merged_df)) "IG")]
  
  
  return(list(interactions = interactions,
              mapped_proteins = final_df))

  
}