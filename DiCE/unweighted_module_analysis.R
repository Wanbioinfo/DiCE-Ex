#' Build reproducible undirected graph
#'
#' Creates an igraph object from an edge list with deterministic ordering.
#' Removes duplicate edges and self-loops, and fixes vertex order for reproducibility.
#'
#' @param edges Data frame with columns `Gene1` and `Gene2`.
#' @return igraph object
#' @noRd
make_unWgraph_repro <- function(edges){
  edges <- edges[order(edges$Gene1, edges$Gene2), ]
  g <- graph_from_data_frame(edges, directed = FALSE)
  if ("simplify" %in% names(formals(simplify))) {
    g <- simplify(g)
  } else {
    g <- igraph::simplify(g)
  }
  
  if (!is.null(V(g)$name)) g <- igraph::permute(g, order(V(g)$name))
  return(g)
}



#' Detect DiCE PPI Modules Using Louvain Community Detection
#'
#' Constructs an unweighted PPI network for the DiCE-selected genes using
#' STRING v12 interactions and identifies community modules via the Louvain
#' algorithm.
#'
#' @param gene_list Character vector of genes to retain in the networks.
#' @param interactions Data frame of interactions with `source`, `target`. 
#' @param algorithm Character. Community detection algorithm to use.
#'   Options are "louvain" (default) or "leiden". Louvain is fast and widely used,
#'   while Leiden provides improved partition quality and guarantees well-connected communities.
#' @param resolution Numeric. Resolution parameter controlling the granularity of detected modules.
#'   Higher values lead to more (smaller) modules, while lower values result in fewer (larger) modules.
#'   Default is 1.
#' @param leiden_itrs Integer. Number of iterations for the Leiden algorithm.
#'   More iterations can improve stability and quality of community detection,
#'   but increase computation time. Default is 3.
#' @param leiden_beta Numeric. Randomness parameter for the Leiden algorithm.
#'   Controls the level of randomness during node reassignment.
#'   Higher values introduce more randomness, potentially escaping local optima.
#'   Default is 0.01.
#' @param seed Integer seed for reproducible community detection. Default is 123.
#' @return A list with:
#'   \itemize{
#'     \item \code{summary_df}: Summary of number of modules and modularity
#'     \item \code{module_stats_df}: Number of nodes and intra-module edges in each module
#'     \item \code{membership_df}: Module assignments and within-module degrees
#'     \item \code{between_module_edges_df}: Number of inter-module edges connecting each pair of modules
#'     \item \code{all_edges_df}: All retained edges with Gene1, Gene2, Gene1_Module, and Gene2_Module
#'     \item \code{edges_by_module}: List of intra-module interaction tables for each module
#'   }
#'
#'@examples
#' \dontrun{
#' # Run module detection on DiCE genes
#' modules <- detect_DiCE_PPI_unweightedModules(
#'     gene_list = my_genes,
#'     interactions = interactions_df,
#'     algorithm = "louvain",
#'     resolution = 0.9,
#'     seed = 123
#' )
#'
#' # View module membership
#' head(modules$membership_df)
#' }
#' @export
detect_PPI_unweightedModules <- function(gene_list = list(), 
                                         interactions = NULL,
                                         algorithm = "louvain",
                                         resolution = 1,
                                         leiden_itrs = 3,
                                         leiden_beta = 0.01,
                                         seed = 123)
{
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  
  if (length(gene_list) == 0){
    stop("Missing gene list!")
  }
  
  interactions <- dplyr::rename(interactions, Gene1 = source, Gene2 = target)
  
  # Keep only interactions where both genes are in the given gene_list
  interactions <- subset(interactions,
                            Gene1 %in% gene_list & Gene2 %in% gene_list)
  
  interactions_unweighted <- interactions[,c("Gene1","Gene2")]
  
  # 1) Build unweighted graph
  # g <- make_unWgraph_repro(interactions_unweighted)
  g <- graph_from_data_frame(interactions_unweighted, directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  if (tolower(algorithm) == "louvain"){
    # 2) Louvain community detection (unweighted)
    cl_obj <- cluster_louvain(graph = g,
                              resolution = resolution)
    memb <- membership(cl_obj)
    
  }else if (tolower(algorithm) == "leiden"){
    # 2) Leiden community detection (unweighted)
    cl_obj <- cluster_leiden(graph = g,
                            objective_function = "modularity",   
                            resolution = resolution,
                            n_iterations = leiden_itrs,
                            beta = leiden_beta)
    memb <- membership(cl_obj)
    
  }else{
    stop("Invalid alogrithm for community detection. Use 'louvain', 'leiden'.")
  }
  
  # Original module labels
  old_module <- as.integer(memb)
  names(old_module) <- names(memb)
  old_module_chr <- as.character(old_module)
  names(old_module_chr) <- names(old_module)
  
  # 3) Annotate edges with module IDs and keep only intra-module edges

  # Temporary edge df with old module names
  edge_df_tmp <- interactions[, c("Gene1", "Gene2"), drop = FALSE]
  edge_df_tmp$Module_Gene1 <- old_module_chr[edge_df_tmp$Gene1]
  edge_df_tmp$Module_Gene2 <- old_module_chr[edge_df_tmp$Gene2]
  
  # Intra-module edges using old modules
  intra_df_tmp <- subset(edge_df_tmp, Module_Gene1 == Module_Gene2)
  intra_df_tmp$Module <- as.character(intra_df_tmp$Module_Gene1)
  
  
  # Nodes per old module
  nodes_per_module_tmp <- data.frame(
    Module = as.character(names(table(old_module_chr))),
    Num_Nodes = as.integer(table(old_module_chr)),
    stringsAsFactors = FALSE
  )
  
  # Edges per old module
  edges_per_module_tmp <- intra_df_tmp %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(Num_Edges = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(Module = as.character(Module))

  
  # Calculate observed / expected ratio
  module_order_df <- dplyr::left_join(
    nodes_per_module_tmp,
    edges_per_module_tmp,
    by = "Module"
  ) %>%
    dplyr::mutate(
      Num_Edges = ifelse(is.na(Num_Edges), 0, Num_Edges),
      Expected_Edges = Num_Nodes * (Num_Nodes - 1) / 2,
      OE_ratio = ifelse(Expected_Edges > 0, Num_Edges / Expected_Edges, 0)
    ) %>%
    dplyr::arrange(
      dplyr::desc(OE_ratio),
      dplyr::desc(Num_Edges),
      dplyr::desc(Num_Nodes)
    ) %>%
    dplyr::mutate(New_Module = paste0("M", dplyr::row_number()))
  
  
  # Old-to-new module map
  module_map <- setNames(module_order_df$New_Module, module_order_df$Module)
  
  # Annotate all edges with NEW module IDs
  edge_df <- interactions[, c("Gene1", "Gene2"), drop = FALSE]
  edge_df$Module_Gene1 <- module_map[old_module_chr[edge_df$Gene1]]
  edge_df$Module_Gene2 <- module_map[old_module_chr[edge_df$Gene2]]
  
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
    Algorithm  = ifelse(tolower(algorithm) == "louvain", "Louvain", "Leiden"),
    Num_Modules  = length(igraph::sizes(cl_obj)),
    Modularity = igraph::modularity(g, memb),
    stringsAsFactors = FALSE
  )
  
  # Membership df using NEW module IDs
  membership_df <- data.frame(
    Gene.Symbol = names(memb),
    Module = module_map[old_module_chr],
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
  module_stats_df <- module_order_df %>%
    dplyr::select(
      Module = New_Module,
      Num_Nodes,
      Num_Edges,
      Expected_Edges,
      OE_ratio
    )
  
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