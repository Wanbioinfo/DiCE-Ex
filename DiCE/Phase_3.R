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
calculate_correlation <- function(class_geneExp, corr_method, corr_mode, corr_pval_cutoff){
  if(corr_mode == "directCorr"){
    corr_matrix <- traditional_corr(class_geneExp, corr_method, corr_pval_cutoff)
    
  }else if (corr_mode == "remove_zerocells"){
    corr_matrix <- remove_Zerocells_corr(class_geneExp, corr_method, corr_pval_cutoff)
    
  }else if (tolower(corr_mode) %in% c("zinbwave", "newwave")){
    corr_matrix <- traditional_corr(class_geneExp, corr_method, corr_pval_cutoff)
    
  }else{
    stop("Invalid correlation mode OR correlation method!")
  }
  
  return(corr_matrix)
}

map_corr_weights <- function(interactions, class_geneExp, corr_matrix){
  colnames(interactions)=c("node1","node2")
  
  vertices <- c(interactions$node1, interactions$node2)
  vertices <- unique(vertices)
  
  # Which network nodes have expression data
  nodes_with_exp <- intersect(vertices,colnames(class_geneExp))
  # List of vertices that do not appear in the expression data
  diff <- setdiff(vertices,nodes_with_exp)
  # Remove those unmatched
  vertices <- vertices[!vertices %in% c(diff)]
  
  # Extract a submatrix from the correlation matrix mat1, limited to only the vertex genes
  vertices <- intersect(vertices, rownames(corr_matrix))
  vertex_corrMat <- corr_matrix[vertices,vertices]
  
  # Melt the correlation matrix into long format
  melt_vertex_corrMat <- stack(as.data.frame(vertex_corrMat))
  melt_vertex_corrMat$ind <- as.character(melt_vertex_corrMat$ind)
  
  ind2 <- rep(row.names(vertex_corrMat), times = ncol(vertex_corrMat))
  melt_vertex_corrMat <- cbind(melt_vertex_corrMat,ind2)
  melt_vertex_corrMat <- melt_vertex_corrMat[order(match(melt_vertex_corrMat[,3],
                                                         interactions[,1])),]
  # Concatenate interacting gene names (node2 + node1)
  concat_pairs <- paste(interactions[,2],interactions[,1])
  
  # Filter matching known STRING interactions
  corr_known_interactions <- melt_vertex_corrMat[paste(melt_vertex_corrMat$ind,
                                                       melt_vertex_corrMat$ind2) %in% concat_pairs,]
  corr_known_interactions <- cbind(corr_known_interactions[3],
                                   corr_known_interactions[2],
                                   corr_known_interactions[1])
  colnames(corr_known_interactions) <- c("source", "target", "weight")
  
  return(corr_known_interactions)
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
                       corr_method, corr_mode, corr_pval_cutoff){
  
  class_geneExp_treatment <- geneExp_to_corr[geneExp_to_corr$class==treatment,]
  class_geneExp_treatment <- class_geneExp_treatment[,-ncol(geneExp_to_corr)]
  
  class_geneExp_control <- geneExp_to_corr[geneExp_to_corr$class==control,]
  class_geneExp_control <- class_geneExp_control[,-ncol(geneExp_to_corr)]
  
  gene_corr_treatment <- calculate_correlation(class_geneExp_treatment, corr_method, corr_mode, corr_pval_cutoff)
  gene_corr_control <- calculate_correlation(class_geneExp_control, corr_method, corr_mode, corr_pval_cutoff)
  
  treatment_ppi_interactions <- map_corr_weights(interactions,
                                               class_geneExp_treatment,
                                               gene_corr_treatment)

  control_ppi_interactions <- map_corr_weights(interactions,
                                             class_geneExp_control,
                                             gene_corr_control)
  

  # Combine interactions edge weights in treatment and control
  merge_interactions_df <- merge(treatment_ppi_interactions, control_ppi_interactions, 
                                 by = c("source", "target"), 
                                 suffixes = c("_treatment", "_control"), all = TRUE)
  
  return(merge_interactions_df)
  
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
#' @param stringDB_confidence Cutoff for confidence score in Stringdb PPI (Default is 400).
#' @return Dataframe of PPIs 
#' @return Dataframe of vertices with the Gene.Symbol and String_ID
#' @noRd
extract_stringdb_PPI <- function(string_protInfo_file, string_ppi_file, stringDB_confidence = 400, phase2_res){
  return(create_stringPPI_fromPhase2(string_protInfo_file, string_ppi_file, stringDB_confidence, phase2_res))
}

#' Run Call create_PPI: Call the function to create PPI from Phase2 genes
#' Helper function (not for users)
#'
#' @param biogrid_ppi_file Filepath of the BioGrid PPI file
#' @param taxonID taxonomy id
#' @param phase2_res Dataframe of Phase 2 results
#' @return Dataframe of PPIs 
#' @return Dataframe of vertices with the Gene.Symbol and String_ID
#' @noRd
extract_biogrid_PPI <- function(biogrid_ppi_file, taxonID, phase2_res){
  return(create_biogridPPI_fromPhase2(biogrid_ppi_file, taxonID, phase2_res))
}
