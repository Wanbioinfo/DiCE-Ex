# ------------------ Phase 3: --------------------------------------
# Making weighted PPI for each phenotype



#' Run calculate_correlation: calculate the correlation 
#' Helper function (not for users)
#'
#' @param class_geneExp Dataframe of gene expression data for specific group of samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#' @param corr_mode Correlation mode ("directCorr" OR "remove_Zerocells" OR "ZINB-WaVE")
#'
#' @return Correlation matrix of the Phase 2 genes
#' @noRd
calculate_correlation <- function(class_geneExp, corr_method, corr_mode){
  if(corr_mode == "directCorr"){
    corr_matrix <- traditional_corr(class_geneExp, corr_method)
    
  }else if (corr_mode == "remove_Zerocells"){
    corr_matrix <- remove_Zerocells_corr(class_geneExp, corr_method)
    
  }else if (corr_mode == "ZINB-WaVE"){
    corr_matrix <- traditional_corr(class_geneExp, corr_method)
    
  }else{
    stop("Invalid correlation mode OR correlation method!")
  }
  
  return(corr_matrix)
}


#' Run Phase 3: Calculate network centralities
#' Helper function (not for users)
#'
#' @param interactions Dataframe of PPI among the Phase 2 genes
#' @param geneExp_to_corr Dataframe of Gene expressions of genes in PPI
#' @param treatment Name of the treatment samples
#' @param control Name of the control samples
#' @param corr_method Correlation method ("pearson" OR "spearman")
#' @param corr_mode Correlation mode ("directCorr" OR "remove_Zerocells" OR "ZINB-WaVE")
#' @param centrality_list Character vector of centrality metrics to compute
#'  
#' @return Dataframe of network centrality values of the genes in the PPI network
#' @return Dataframe of edges with the corresponding weights (correlation coefficient) 
#' @noRd
run_phase3 <- function(interactions, geneExp_to_corr, treatment, control, 
                       corr_method, corr_mode, centrality_list){
  
  class_geneExp_treatment <- geneExp_to_corr[geneExp_to_corr$class==treatment,]
  class_geneExp_treatment <- class_geneExp_treatment[,-ncol(geneExp_to_corr)]
  
  class_geneExp_control <- geneExp_to_corr[geneExp_to_corr$class==control,]
  class_geneExp_control <- class_geneExp_control[,-ncol(geneExp_to_corr)]
  
  gene_corr_treatment <- calculate_correlation(class_geneExp_treatment, corr_method, corr_mode)
  gene_corr_control <- calculate_correlation(class_geneExp_control, corr_method, corr_mode)
  
  treatment_ppi_output <- PPI_network_analysis(interactions,
                                               class_geneExp_treatment,
                                               gene_corr_treatment,
                                               centrality_list)
  centralities_treatment <- treatment_ppi_output$centralities
  interactions_treatment <- treatment_ppi_output$corr_known_interactions
  
  # add suffix
  centralities_treatment <- as.data.frame(centralities_treatment)
  suffix <- "_treatment" 
  colnames(centralities_treatment)[-1] <- paste0(colnames(centralities_treatment)[-1], suffix)
  
  
  control_ppi_output <- PPI_network_analysis(interactions,
                                             class_geneExp_control,
                                             gene_corr_control,
                                             centrality_list)
  
  centralities_control <- control_ppi_output$centralities
  interactions_control <- control_ppi_output$corr_known_interactions
  
  # add suffix
  centralities_control <- as.data.frame(centralities_control)
  suffix <- "_control" 
  colnames(centralities_control)[-1] <- paste0(colnames(centralities_control)[-1], suffix)
  
  
  # Merge treatment and control centralities
  merge_centralities_df <- merge(centralities_treatment, centralities_control, by = "Gene.Symbol")

  # Combine interactions edge weights in treatment and control
  merge_interactions_df <- merge(interactions_treatment, interactions_control, 
                                 by = c("source", "target"), 
                                 suffixes = c("_treatment", "_control"), all = TRUE)
  
  return(list(phase3_centralities_df = merge_centralities_df,
              interactions_withWeights_df = merge_interactions_df))
  
}

#' Run Get network gene expressions: Filter expression data to only include genes that are in the PPI.
#' Helper function (not for users)
#'
#' @param geneExp_data Dataframe of gene expressions of all genes
#' @param mapped_proteins Dataframe of genes in the PPI network
#'
#' @return Dataframe of filtered gene expressions
#' @noRd
getNetworkGene_expression <- function(geneExp_data,mapped_proteins){
  class <- geneExp_data[,ncol(geneExp_data)]
  filtered_geneExp_data <- geneExp_data[,colnames(geneExp_data) %in% mapped_proteins$Gene.Symbol]
  filtered_geneExp_data <- as.data.frame(filtered_geneExp_data)
  
  filtered_geneExp_data <- cbind(filtered_geneExp_data, class)
  
  #table(duplicated(colnames(filtered_geneExp_data)))
  
  return(filtered_geneExp_data)
}

#' Run Call create_PPI: Call the function to create PPI from Phase2 genes
#' Helper function (not for users)
#'
#' @param string_protInfo_file Filepath of the StringDB protein information file
#' @param string_ppi_file Filepath of the StringDB PPI file
#' @param phase2_res Dataframe of Phase 2 results
#'
#' @return Dataframe of PPIs 
#' @return Dataframe of vertices with the Gene.Symbol and String_ID
#' @noRd
extract_PPI <- function(string_protInfo_file, string_ppi_file, phase2_res){
  return(create_PPI_fromPhase2(string_protInfo_file, string_ppi_file, phase2_res))
}
