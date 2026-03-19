#' Detect DiCE PPI Modules Using Louvain Community Detection
#'
#' Constructs an unweighted PPI network for the DiCE-selected genes using
#' STRING v12 interactions and identifies community modules via the Louvain
#' algorithm.
#'
#' @param gene_list Character vector of genes to retain in the networks.
#' @param interactions_df Data frame of interactions with `source`, `target`. 
#' @param seed Integer seed for reproducible community detection. Default is 123.
#' @return A list with:
#'   \itemize{
#'     \item \code{summary_df}: Summary of number of modules and modularity
#'     \item \code{module_stats_df}: Number of nodes and intra-module edges in each module
#'     \item \code{membership_df}: Module assignments and within-module degrees
#'     \item \code{between_module_edges_df}: Number of inter-module edges connecting each pair of modules
#'     \item \code{all_edges_df}: All retained edges with Gene1, Gene2, combined_score, Gene1_Module, and Gene2_Module
#'     \item \code{edges_by_module}: List of intra-module interaction tables for each module
#'   }
#'
#'@examples
#' \dontrun{
#' # Run module detection on DiCE genes
#' modules <- detect_DiCE_PPI_unweightedModules(
#'     gene_list = my_genes,
#'     species = "human",
#'     seed = 123
#' )
#'
#' # View module membership
#' head(modules$membership_df)
#' }
#' @export
detect_PPI_unweightedModules <- function(gene_list = list(), 
                                         interactions_df = NULL,
                                         seed = 123)
{
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  
  if (length(gene_list) == 0){
    stop("Missing gene list!")
  }
  
  interactions_df <- dplyr::rename(interactions_df, Gene1 = source, Gene2 = target)
  
  # Keep only interactions where both genes are in the given gene_list
  interactions_df <- subset(interactions_df,
                              Gene1 %in% gene_list & Gene2 %in% gene_list)
  
  interactions_unweighted <- interactions_df[,c("Gene1","Gene2")]
  
  # 1) Build unweighted graph
  g <- graph_from_data_frame(interactions_unweighted, directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  # 2) Louvain community detection (unweighted)
  cl_louvain <- cluster_louvain(g)
  memb <- membership(cl_louvain)
  
  # 3) Annotate edges with module IDs and keep only intra-module edges

  # Annotate all edges with module IDs
  edge_df <- interactions_df[, c("Gene1", "Gene2"), drop = FALSE]
  edge_df$Module_Gene1 <- paste0("M", as.integer(memb[edge_df$Gene1]))
  edge_df$Module_Gene2 <- paste0("M", as.integer(memb[edge_df$Gene2]))
  
  # Intra-module edges only
  intra_df <- subset(edge_df, Module_Gene1 == Module_Gene2)
  intra_df$Module <- intra_df$Module_Gene1
  
  # Split intra-module interactions per module
  edges_by_module <- split(
    intra_df[, c("Gene1", "Gene2"), drop = FALSE],
    intra_df$Module
  )
  
  # Summary df
  summary_df <- data.frame(
    Algorithm  = "Louvain",
    Num_Modules  = length(igraph::sizes(cl_louvain)),
    Modularity = igraph::modularity(g, memb),
    stringsAsFactors = FALSE
  )
  
  # Membership df
  membership_df <- data.frame(
    Gene.Symbol   = names(memb),
    Module = paste0("M", as.integer(memb)),
    stringsAsFactors = FALSE
  )
  
  # Degree inside module only
  membership_df$Degree_inModule <- vapply(seq_along(memb), function(i) {
    node <- names(memb)[i]
    mod  <- memb[[i]]
    vs   <- names(memb[memb == mod])
    subg <- igraph::induced_subgraph(g, vids = vs)
    as.integer(igraph::degree(subg)[node])
  }, integer(1))
  
  # nodes and edges per module
  nodes_per_module <- membership_df %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(Num_Nodes = dplyr::n(), .groups = "drop")
  
  edges_per_module <- intra_df %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(Num_Edges = dplyr::n(), .groups = "drop")
  
  module_stats_df <- dplyr::left_join(nodes_per_module, edges_per_module, by = "Module")
  module_stats_df$Num_Edges[is.na(module_stats_df$Num_Edges)] <- 0
  
  # edges among modules
  inter_df <- subset(edge_df, Module_Gene1 != Module_Gene2)
  
  if (nrow(inter_df) > 0) {
    inter_df$Module_A <- pmin(inter_df$Module_Gene1, inter_df$Module_Gene2)
    inter_df$Module_B <- pmax(inter_df$Module_Gene1, inter_df$Module_Gene2)
    
    between_module_edges_df <- inter_df %>%
      dplyr::group_by(Module_A, Module_B) %>%
      dplyr::summarise(Num_Edges = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(Module_A, Module_B)
  } else {
    between_module_edges_df <- data.frame(
      Module_A = character(0),
      Module_B = character(0),
      `Num_Edges (between modules)` = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    summary_df = summary_df,
    module_stats_df = module_stats_df,
    membership_df = membership_df,
    between_module_edges_df = between_module_edges_df,
    all_edges_df = edge_df,
    edges_by_module = edges_by_module
  ))
}


