#' Detect and compare weighted PPI network modules across two conditions
#'
#' Builds condition-specific weighted PPI networks (treatment vs control) from a
#' shared interaction table, detects Louvain communities, summarizes modules, and
#' compares module overlap and node-level module roles between conditions.
#'
#' @param gene_list Character vector of genes to retain in the networks.
#' @param interactions_df Data frame of interactions with `source`, `target`, and
#' condition-specific weight columns named `weight_<condition>` (e.g., `weight_Tumor`).
#' @param treatment Character string naming the treatment condition (used to select
#' `weight_<treatment>`. E.g., "Tumor").
#' @param control Character string naming the control condition (used to select
#' `weight_<control>`. E.g., "Normal").
#' @param louvain_resolution  Resolution parameter for computing modularity in Louvain algorithm (default 1).
#' @param seed Random seed for reproducibility (default 123).
#'
#' @return A list with:
#' \describe{
#'   \item{network_modules}{Data frame summarizing community detection per condition (e.g., n_modules, modularity).}
#'   \item{module_stats}{Module-level table (module size, edge count, density, mean weight) for each condition.}
#'   \item{cmp_jaccard}{Matrix/data frame of Jaccard similarity between treatment and control modules.}
#'   \item{cmp_count}{Matrix/data frame of gene overlap counts between treatment and control modules.}
#'   \item{nodes_stats}{Node-level table comparing module/role statistics across conditions (joined by gene).}
#' }
#'
#' @examples
#' \dontrun{
#' res <- detect_PPI_weightedModules(
#'   gene_list = my_genes,
#'   interactions_df = intr_df,
#'   treatment = "Tumor",
#'   control = "Normal",
#'   louvain_resolution = 0.9,
#'   seed = 123
#' )
#' 
#' # View network stats of the nodes in each module
#' head(res$nodes_stats)
#' 
#' }
#'
#' @export
detect_PPI_weightedModules <- function(gene_list = list(),
                                       interactions_df = NULL,
                                       treatment = NULL,
                                       control = NULL,
                                       louvain_resolution = 1,
                                       seed = 123)
{
  
  if (length(gene_list) == 0){
    stop("Missing gene list!")
  }
  
  if (is.null(interactions_df)){
    stop("Missing interactions data!")
  }
  
  interactions_df <- dplyr::rename(interactions_df, Gene1 = source, Gene2 = target)
  
  # Keep only interactions where both genes are in the given gene_list
  keep_interactions <- subset(interactions_df,
                         Gene1 %in% gene_list & Gene2 %in% gene_list)
  
  # create weighted graphs separately for treatment and control data
  treatment_intrs <- filter_weightedEdges(keep_interactions, treatment)
  control_intrs <- filter_weightedEdges(keep_interactions, control)
  
  # Detect network communities for treatment and control networks
  treatment_comDetection <- detect_communities(treatment_intrs,louvain_resolution)
  treatment_graph <- treatment_comDetection$graph
  treatment_communities <- treatment_comDetection$communities
  treatment_membs <- treatment_comDetection$membership
  
  control_comDetection <- detect_communities(control_intrs,louvain_resolution)
  control_graph <- control_comDetection$graph
  control_communities <- control_comDetection$communities
  control_membs <- control_comDetection$membership

  # Get treatment and control community details
  treatment_comDetails <- extract_communityDetails(treatment, treatment_intrs,treatment_graph,
                                                   treatment_membs,treatment_communities)
  treatment_summary_df <- treatment_comDetails$summary_df
  treatment_memb_df <- treatment_comDetails$membership_df
  treatment_edges_with_module <- treatment_comDetails$edges_with_module
  
  control_comDetails <- extract_communityDetails(control, control_intrs,control_graph,
                                                 control_membs,control_communities)
  control_summary_df <- control_comDetails$summary_df
  control_memb_df <- control_comDetails$membership_df
  control_edges_with_module <- control_comDetails$edges_with_module
  
  treatment_summary_df$Condition <- treatment
  control_summary_df$Condition   <- control
  
  combined_summary_df <- dplyr::bind_rows(treatment_summary_df, control_summary_df)
  
  # Module summaries
  treatment_modSummary <- summarize_modules(treatment_comDetection, treatment)
  control_modSummary <- summarize_modules(control_comDetection, control)
  
  combined_modSummary <- dplyr::bind_rows(treatment_modSummary, control_modSummary)
  
  # Analysis 1 - Calculate jaccard score between modules
  cmp1 <- compare_ConditionModules(treatment_membs, control_membs)
  
  cmp_jaccard <- as.data.frame(cmp1$jaccard)
  rownames(cmp_jaccard) <- sub("^Treatment", treatment, rownames(cmp_jaccard))
  colnames(cmp_jaccard) <- sub("^Control", control, colnames(cmp_jaccard))
  
  cmp_count <- as.data.frame(cmp1$count)
  rownames(cmp_count) <- sub("^Treatment", treatment, rownames(cmp_count))
  colnames(cmp_count) <- sub("^Control", control, colnames(cmp_count))
  
  
  # Analysis 2 - Node features in each module
  nodes_treatment <- node_moduleInfo(treatment_comDetection, treatment, weight_attr = "weight")
  nodes_control <- node_moduleInfo(control_comDetection, control, weight_attr = "weight")
  
  nodes_combined <- full_join(nodes_treatment, nodes_control, by = "Gene.Symbol", 
                           suffix = c(paste0("_",treatment), paste0("_",control)))
  
  
  # Analysis 3 - get all interactions with modules information and edge weights
  edges_with_modules <- full_join(treatment_edges_with_module, control_edges_with_module, by = c("Gene1", "Gene2"))
  
  edges_with_modules <- edges_with_modules[, c(
    "Gene1",paste0("Module_Gene1_",treatment), paste0("Module_Gene1_",control),
    "Gene2",paste0("Module_Gene2_",treatment),paste0("Module_Gene2_",control),
    paste0("weight_",treatment),paste0("weight_",control)
    
  )]
  
  # Analysis 4 - condition specific intermodule connectivity
  treatment_interMod <- inter_module_connectivity(edges_with_modules, treatment)
  control_interMod <- inter_module_connectivity(edges_with_modules, control)
  
  return(list(network_modules = combined_summary_df,
              module_stats = combined_modSummary,
              cmp_jaccard = cmp_jaccard,
              cmp_count = cmp_count,
              nodes_stats = nodes_combined,
              edges_with_modules = edges_with_modules,
              treatment_interMod = treatment_interMod,
              control_interMod = control_interMod
              ))

  
}