#' Build reproducible undirected graph
#'
#' Creates an igraph object from an edge list with deterministic ordering.
#' Removes duplicate edges and self-loops, and fixes vertex order for reproducibility.
#'
#' @param edges Data frame with columns `Gene1`, `Gene2` and `weight`.
#' @return igraph object
#' @noRd
make_graph_repro <- function(edges){
  edges <- edges[order(edges$Gene1, edges$Gene2), ]
  
  # Build graph with weights (assumes column name = "weight")
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  
  # Ensure weight attribute exists
  if (!"weight" %in% colnames(edges)) {
    stop("Input edges must contain a 'weight' column.")
  }
  
  E(g)$weight <- edges$weight
  
  # Simplify: combine duplicate edges by taking mean weights
  if ("simplify" %in% names(formals(simplify))) {
    g <- simplify(g, edge.attr.comb = list(weight = "mean"))
  } else {
    g <- igraph::simplify(g, edge.attr.comb = list(weight = "mean"))
  }
  
  # Reproducible vertex order
  if (!is.null(V(g)$name)) {
    g <- igraph::permute(g, order(V(g)$name))
  }
  
  return(g)
}


#' Extract condition-specific weighted edges
#'
#' Internal helper to select the edge weight column for a given condition
#' (`weight_<condition>`) and return a standardized edge list.
#'
#' @param all_edges Data frame with `Gene1`, `Gene2`, and weight columns.
#' @param condition Condition name used to select the weight column.
#'
#' @return Data frame with columns `Gene1`, `Gene2`, and `weight`.
#' @noRd
filter_weightedEdges <- function(all_edges, condition) {

  wcol <- paste0("weight_", condition)
  
  stopifnot(all(c("Gene1","Gene2", wcol) %in% names(all_edges)))
  
  weights <- all_edges[[wcol]]
  
  all_edges <- all_edges %>% transmute(Gene1, Gene2, weight = as.numeric(weights))
  
  return(all_edges)
}

#' Detect PPI network communities using Louvain
#'
#' Internal helper to construct a weighted graph from interaction data and
#' identify network modules using the Louvain community detection algorithm.
#'
#' @param interactions Data frame with `Gene1`, `Gene2`, and `weight`.
#' @param louvain_resolution  Resolution parameter for computing modularity in Louvain algorithm (default 1).
#' @param seed Random seed for reproducibility.
#'
#' @return List containing the igraph object, detected communities, and node membership.
#' @noRd
detect_communities <- function(interactions,
                               algorithm,
                               resolution,
                               leiden_itrs,
                               leiden_beta,
                               seed){
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  
  # convert the weight for the absolute value
  interactions$weight <- abs(as.numeric(interactions$weight))
  
  # Build weighted graph
  # graph <- make_graph_repro(interactions)

  graph <- graph_from_data_frame(interactions, directed = FALSE)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
  
  # Community detection
  if (tolower(algorithm) == "louvain"){
    # 2) Louvain community detection (unweighted)
    cl_obj <- cluster_louvain(graph = graph,
                              weights = E(graph)$weight,
                              resolution = resolution)
    memb <- membership(cl_obj)
    
  }else if (tolower(algorithm) == "leiden"){
    # 2) Leiden community detection (unweighted)
    cl_obj <- cluster_leiden(graph = graph,
                             weights = E(graph)$weight,
                             objective_function = "modularity",   
                             resolution = resolution,
                             n_iterations = leiden_itrs,
                             beta = leiden_beta)
    memb <- membership(cl_obj)
    
  }else{
    stop("Invalid alogrithm for community detection. Use 'louvain', 'leiden'.")
  }
  
  
  return(list(graph = graph,
              communities = cl_obj,
              membership = memb))
  
}

#' Extract community module details
#'
#' Internal helper to summarize detected communities, compute node
#' membership and within-module degree, and split intra-module edges.
#'
#' @param condition Character string specifying the condition name.
#' @param interactions Data frame with `Gene1`, `Gene2`, and `weight`.
#' @param graph igraph object of the network.
#' @param membership Named vector of node module assignments.
#' @param communities Community object returned by igraph.
#'
#' @return List containing module summary, node membership table,
#' and intra-module edges grouped by module.
#' @noRd
extract_communityDetails <- function(condition, interactions,graph, membership, communities){
  # Annotate edges with module IDs and keep only intra-module edges
  edge_cols <- intersect(c("Gene1","Gene2","weight"), colnames(interactions))
  edge_df <- interactions[, edge_cols, drop = FALSE]
  edge_df$Module_Gene1 <- membership[edge_df$Gene1]
  edge_df$Module_Gene2 <- membership[edge_df$Gene2]
  
  edges_with_modules <- edge_df[, c("Gene1","Gene2","weight",
                                    "Module_Gene1","Module_Gene2")]
  
  edges_with_modules <- edges_with_modules %>%
    dplyr::mutate(
      Module_Gene1 = paste0(condition, "_M", Module_Gene1),
      Module_Gene2 = paste0(condition, "_M", Module_Gene2)
    ) %>%
    dplyr::rename(
      !!paste0("weight_", condition) := weight,
      !!paste0("Module_Gene1_", condition) := Module_Gene1,
      !!paste0("Module_Gene2_", condition) := Module_Gene2
    )

  intra_df <- subset(edge_df, Module_Gene1 == Module_Gene2)
  intra_df$Module <- intra_df$Module_Gene1
  intra_df$Module_Gene1 <- NULL
  intra_df$Module_Gene2 <- NULL
  
  # Split intra-module interactions per module
  edges_by_module <- split(intra_df[, edge_cols, drop = FALSE], intra_df$Module)
  
  # Prepare summary + membership
  summary_df <- data.frame(
    Algorithm  = "Louvain",
    n_modules  = length(sizes(communities)),
    Modularity = modularity(graph, membership),
    stringsAsFactors = FALSE
  )
  
  membership_df <- data.frame(
    Gene   = names(membership),
    Module = as.integer(membership),
    stringsAsFactors = FALSE
  )
  
  # Degree inside module only
  membership_df$Degree_inModule <- vapply(seq_along(membership), function(i) {
    node <- names(membership)[i]
    mod  <- membership[[i]]
    vs   <- names(membership[membership == mod])
    subg <- induced_subgraph(graph, vids = vs)
    as.integer(degree(subg)[node])
  }, integer(1))
  
  return(list(summary_df = summary_df,
              membership_df = membership_df,
              edges_with_modules = edges_with_modules))
}

#' Calculate module density statistics
#'
#' Internal helper to compute module size, number of edges, unweighted
#' density, and mean edge weight for a given module subgraph.
#'
#' @param g igraph object representing the network.
#' @param module_genes Character vector of genes in the module.
#'
#' @return Tibble with module size, edge count, density, and mean weight.
#' @noRd
module_stats <- function(g, module_genes) {
  sg <- induced_subgraph(g, vids = intersect(module_genes, V(g)$name))
  vertices  <- vcount(sg)
  edges <- ecount(sg)
  if (vertices < 2) {
    return(tibble(Module_Size=vertices, 
                  Num_Edges=edges, 
                  Density_Unweighted = NA_real_, 
                  Mean_Weight = NA_real_))
  }
  density_unw <- (2*edges) / (vertices*(vertices-1))
  w <- E(sg)$weight
  tibble(Module_Size=vertices, 
         Num_Edges=edges, 
         Density_Unweighted = density_unw, 
         Mean_Weight = if (length(w)) mean(w) else NA_real_)
}

#' Summarize modules for a condition
#'
#' Internal helper to compute module size and density statistics for
#' detected communities and label modules with a condition-specific prefix.
#'
#' @param res List containing graph and community detection results.
#' @param condition Character string specifying the condition name.
#'
#' @return Tibble summarizing modules with size and density metrics.
#' @noRd
summarize_modules <- function(res, condition) {
  g    <- res$graph
  memb <- membership(res$communities)
  mods <- split(names(memb), memb)
  
  # label modules with condition-specific prefix
  prefix <- paste0(condition, "_M")
  module_names <- paste0(prefix, names(mods))
  
  tibble::tibble(
    Module = module_names
  ) %>%
    dplyr::mutate(tmp = purrr::map(mods, ~ module_stats(g, .x))) %>%
    tidyr::unnest_wider(tmp) %>%
    dplyr::mutate(Condition = condition)
}

#' Compute Jaccard similarity between modules
#'
#' Internal helper to calculate gene overlap count and Jaccard similarity
#' between treatment and control modules.
#'
#' @param mods_t List of treatment modules (gene vectors).
#' @param mods_c List of control modules (gene vectors).
#'
#' @return List containing Jaccard similarity and overlap count matrices.
#' @noRd
jaccard_modules <- function(mods_t, mods_c) {
  Mt <- length(mods_t); Mn <- length(mods_c)
  J  <- matrix(0, nrow = Mt, ncol = Mn,
               dimnames = list(paste0("Treatment_M", names(mods_t)),
                               paste0("Control_M", names(mods_c))))
  
  C  <- matrix(0, nrow = Mt, ncol = Mn,
                   dimnames = list(paste0("Treatment_M", names(mods_t)),
                                   paste0("Control_M", names(mods_c))))
  
  
  for (i in seq_len(Mt)) for (j in seq_len(Mn)) {
    a <- mods_t[[i]]; b <- mods_c[[j]]
    C[i,j] <- length(intersect(a,b))
    J[i,j] <- length(intersect(a,b)) / length(union(a,b))
  }
  
  return(list(jaccard = J,
              count = C))
}

#' Compare treatment and control modules
#'
#' Internal helper to compute gene overlap and Jaccard similarity
#' between treatment and control network modules.
#'
#' @param treatment_membs Named vector of treatment module memberships.
#' @param control_membs Named vector of control module memberships.
#'
#' @return List containing Jaccard similarity and overlap count matrices.
#' @noRd
compare_ConditionModules <- function(treatment_membs, control_membs) {
  # Build module membership lists
  mods_t <- split(names(treatment_membs), treatment_membs)   # Treatment modules
  mods_c <- split(names(control_membs),   control_membs)     # Control modules
  
  # Jaccard matrix (rows: Treatment, cols: Control)
  cmp <- jaccard_modules(mods_t, mods_c)
  jaccard <- cmp$jaccard
  count <- cmp$count
  
  return(list(jaccard = jaccard,
              count = count))
}

#' Calculate node role statistics in a weighted network
#'
#' Internal helper to compute each node’s within-module z-score (from weighted
#' within-module degree) and weighted participation coefficient across modules.
#'
#' @param g igraph object representing the network.
#' @param comm Community object returned by igraph.
#' @param condition Character string used to prefix module labels.
#' @param weight_attr Edge attribute name storing weights (default `"weight"`).
#'
#' @return Tibble with `Gene.Symbol`, `Module`, `k_in`, `k_out`, `k_total`, `z`, and `P`.
#' @noRd
calculate_node_stats <- function(g, comm, condition, weight_attr = "weight") {
  
  memb <- membership(comm)
  
  # assign module labels
  V(g)$module <- paste0(condition,"_M", memb[V(g)])

  
  # weighted adjacency (sparse)
  adj <- as_adjacency_matrix(g, sparse = TRUE, attr = weight_attr)
  
  # total weighted degree (strength)
  k_total <- Matrix::rowSums(adj)
  
  # module bookkeeping
  mod_levels <- unique(V(g)$module)
  mod_id <- match(V(g)$module, mod_levels)
  
  # node x module indicator (sparse)
  M <- sparseMatrix(
    i = seq_len(vcount(g)),
    j = mod_id,
    x = 1,
    dims = c(vcount(g), length(mod_levels))
  )
  
  # k_by_mod[i, m] = total edge weight from node i to module m
  k_by_mod <- adj %*% M   # (n x n) %*% (n x M) => (n x M)
  
  # within-module weighted degree
  k_in <- as.numeric(k_by_mod[cbind(seq_len(vcount(g)), mod_id)])
  
  # outside-module weighted degree
  k_out <- as.numeric(k_total - k_in)
  
  # within-module degree z-score (computed from weighted k_in)
  df0 <- tibble(
    Gene.Symbol    = V(g)$name,
    Module  = V(g)$module,
    k_in    = k_in,
    k_out   = k_out,
    k_total = as.numeric(k_total)
  )
  
  df_z <- df0 %>%
    group_by(Module) %>%
    mutate(k_in_mu = mean(k_in),
           k_in_sd = sd(k_in)) %>%
    ungroup() %>%
    mutate(z = ifelse(k_in_sd > 0, (k_in - k_in_mu) / k_in_sd, 0)) %>%
    dplyr::select(-k_in_mu, -k_in_sd)
  
  # weighted participation coefficient
  frac_sq <- Matrix::rowSums((k_by_mod / pmax(1e-12, k_total))^2)
  P <- as.numeric(1 - frac_sq)
  
  df_z %>% mutate(P = P)
}

#' Extract node module role information
#'
#' Internal helper to compute node-level module statistics and format them
#' into a table with within-module z-score, participation coefficient,
#' and intra/inter-module degree.
#'
#' @param res List containing graph and community detection results.
#' @param condition Character string specifying the condition name.
#' @param weight_attr Edge attribute storing weights (default `"weight"`).
#'
#' @return Tibble with node module role statistics.
#' @noRd
node_moduleInfo <- function(res, condition, weight_attr = "weight") {
  calculate_node_stats(res$graph, res$communities, condition, weight_attr = weight_attr) %>%
    mutate(Condition = condition) %>%
    dplyr::rename(
      Within_Module_Z = z,
      Participation_Coefficient = P,
      Intra_Module_Degree = k_in,
      Inter_Module_Degree = k_out
    ) %>%
    dplyr::select(Condition, Gene.Symbol, Module,
           Within_Module_Z,
           Participation_Coefficient,
           Intra_Module_Degree,
           Inter_Module_Degree)
}

#' Summarize inter-module connectivity within a condition
#'
#' @param all_edges A data frame containing all the edge weights and modules information
#' @param condition Character string specifying the condition name.
#'
#' @return A data frame with one row per inter-module pair
#' @noRd
inter_module_connectivity <- function(edges_with_modules, condition){
  
  mod1_col   <- paste0("Module_Gene1_", condition)
  mod2_col   <- paste0("Module_Gene2_", condition)
  weight_col <- paste0("weight_", condition)
  
  inter_module <- edges_with_modules %>%
    dplyr::filter(.data[[mod1_col]] != .data[[mod2_col]])
  
  inter_module_summary <- inter_module %>%
    dplyr::mutate(
      Module_A = pmin(.data[[mod1_col]], .data[[mod2_col]]),
      Module_B = pmax(.data[[mod1_col]], .data[[mod2_col]])
    ) %>%
    dplyr::group_by(Module_A, Module_B) %>%
    dplyr::summarise(
      edge_count  = dplyr::n(),
      weight_sum  = sum(.data[[weight_col]], na.rm = TRUE),
      weight_mean = mean(.data[[weight_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  return(inter_module_summary)
}


#' Rename modules by observed/expected edge ratio
#'
#' Reorders modules based on intra-module observed-to-expected edge ratio
#' and assigns new module labels (M1 = highest OE).
#'
#' @param module_stats_df Data frame with module size and edge counts.
#' @param condition Condition name prefix.
#' @return List with renamed module stats and old-to-new module map.
#' @noRd
rename_modules_by_OE <- function(module_stats_df, condition) {
  
  module_stats_df <- module_stats_df %>%
    dplyr::mutate(
      Expected_Edges = Module_Size * (Module_Size - 1) / 2,
      OE_ratio = ifelse(Expected_Edges > 0, Num_Edges / Expected_Edges, 0)
    ) %>%
    dplyr::arrange(
      dplyr::desc(OE_ratio),
      dplyr::desc(Num_Edges),
      dplyr::desc(Module_Size)
    ) %>%
    dplyr::mutate(
      Old_Module = Module,
      New_Module = paste0(condition, "_M", dplyr::row_number())
    )
  
  module_map <- setNames(module_stats_df$New_Module, module_stats_df$Old_Module)
  
  list(
    module_stats_df = module_stats_df %>%
      dplyr::select(
        Module = New_Module,
        Condition,
        Module_Size,
        Num_Edges,
        Expected_Edges,
        OE_ratio,
        Density_Unweighted,
        Mean_Weight
      ),
    module_map = module_map
  )
}