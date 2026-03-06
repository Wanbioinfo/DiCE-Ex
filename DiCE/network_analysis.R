# ------------------------- Network analysis for Phase 3 -------------------------------
#' Run PPi network analysis : Calculate betweenness and eigen vector centralities
#' Helper function (not for users)
#'
#' @param interactions Dataframe of PPIs
#' @param class_geneExp Dataframe of gene expression for a specific class of samples
#' @param corr_matrix Matrix of absolute correlation coefficients between genes
#' @param centrality_list Character vector of centrality metrics to compute.
#'
#' @return Dataframe of network centrality values of the genes in the PPI network 
#' @return Dataframe of edges with the corresponding weights (correlation coefficient) 
#' @noRd
PPI_network_analysis <- function(interactions,class_geneExp,corr_matrix,centrality_list){
  
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
  
  # Calculate network centralities
  cent_value_list <- list()
  
  for (centrality in centrality_list){
    centrality = tolower(centrality)
    
    if ((centrality == "betweenness") | (centrality == "betweeness")){
      cent_df <- calc_betweenness(corr_known_interactions)
      
    }else if ((centrality == "eigen vector") | 
              (centrality == "eig") | 
              (centrality == "eigenvector")){
      cent_df <- calc_eigenVC(corr_known_interactions)
      
    }else if (centrality == "authority"){
      cent_df <- calc_authority(corr_known_interactions)
      
    }else if (centrality == "strength"){
      cent_df <- calc_strength(corr_known_interactions)
      
    }else if (centrality == "closeness"){
      cent_df <- calc_closeness(corr_known_interactions)
      
    }else if (centrality == "pagerank"){
      cent_df <- calc_pagerank(corr_known_interactions)
      
    }else if (centrality == "harmonic"){
      cent_df <- calc_harmonic(corr_known_interactions)
      
    }else {
      stop("Invalid network centrality measure or it is not supported by DiCE. 
           DiCE supports for betweeness/eigen vector/authority/strength/closeness/pagerank/harmonic !")
    }
    
    # Convert row names to a 'Gene.Symbol' column
    cent_df <- cent_df %>%
      tibble::rownames_to_column(var = "Gene.Symbol")
    
    cent_value_list[[centrality]] <- cent_df
  }
  
  # Combine all by 'Gene' column
  combined_centralities <- purrr::reduce(cent_value_list, full_join, by = "Gene.Symbol")

  combined_centralities <- combined_centralities %>% relocate(Gene.Symbol)
  
  
  return(list(centralities=combined_centralities,
              corr_known_interactions = corr_known_interactions))
}

