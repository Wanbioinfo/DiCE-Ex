#' Detect DiCE PPI Modules Using Louvain Community Detection
#'
#' Constructs an unweighted PPI network for the DiCE-selected genes using
#' STRING v12 interactions and identifies community modules via the Louvain
#' algorithm.
#'
#' @param gene_list Character vector of genes to retain in the networks.
#' @param species Character string indicating the organism
#'   (\code{"human"} or \code{"mouse"}), used to load STRING reference files. Default = "human".
#' @param seed Integer seed for reproducible community detection. Default is 123.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{summary_df}: Summary of number of modules and modularity
#'     \item \code{membership_df}: Module assignments and within-module degrees
#'     \item \code{edges_by_module}: List of intra-module interaction tables
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
                                         species = "human",
                                         seed = 123)
{
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  
  if (length(gene_list) == 0){
    stop("Missing gene list!")
  }
  
  # String db downloaded files
  if(species == "human"){
    string_protInfo_file <- "extdata/stringDB_v12/human/9606.protein.info.v12.0.txt"
    string_ppi_file <- "extdata/stringDB_v12/human/9606.protein.links.v12.0.txt.gz"
    #    string_protInfo_file <- system.file("extdata/stringDB_v12/human/9606.protein.info.v12.0.txt", package = "DiCE")
    #    string_ppi_file <- system.file("extdata/stringDB_v12/human/9606.protein.links.v12.0.txt", package = "DiCE")
    taxonID <- 9606
  }else if(species == "mouse"){
    string_protInfo_file <- "extdata/stringDB_v12/mouse/10090.protein.info.v12.0.txt"
    string_ppi_file <- "extdata/stringDB_v12/mouse/10090.protein.links.v12.0.txt.gz"
    #    string_protInfo_file <- system.file("extdata/stringDB_v12/mouse/10090.protein.info.v12.0.txt", package = "DiCE")
    #    string_ppi_file <- system.file("extdata/stringDB_v12/mouse/10090.protein.links.v12.0.txt", package = "DiCE")
    taxonID = 10090
  }
  
  
  # Read Protein information data
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
  mapped_proteins <- mapped_proteins %>% filter(Gene.Symbol %in% gene_list)
  
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
  ppi_df <- ppi_df %>% filter(combined_score >= 400)
  
  # Map gene symbols -> STRING IDs
  string_ids <- mapped_proteins$STRING_id[mapped_proteins$Gene.Symbol %in% gene_list]
  
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
  
  interactions_unweighted <- interactions[,c("Gene1","Gene2")]
  
  # 1) Build unweighted graph
  g <- graph_from_data_frame(interactions_unweighted, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  # 2) Louvain community detection (unweighted)
  cl_louvain <- cluster_louvain(g)
  memb <- membership(cl_louvain)
  
  # 3) Annotate edges with module IDs and keep only intra-module edges
  edge_cols <- intersect(c("Gene1","Gene2","combined_score"), colnames(interactions))
  edge_df <- interactions[, edge_cols, drop = FALSE]
  edge_df$Module_Gene1 <- memb[edge_df$Gene1]
  edge_df$Module_Gene2 <- memb[edge_df$Gene2]
  
  intra_df <- subset(edge_df, Module_Gene1 == Module_Gene2)
  intra_df$Module <- intra_df$Module_Gene1
  intra_df$Module_Gene1 <- NULL
  intra_df$Module_Gene2 <- NULL
  
  # 4) Split intra-module interactions per module
  edges_by_module <- split(intra_df[, edge_cols, drop = FALSE], intra_df$Module)
  
  # 5) Prepare summary + membership
  summary_df <- data.frame(
    Algorithm  = "Louvain",
    n_modules  = length(sizes(cl_louvain)),
    Modularity = modularity(g, memb),
    stringsAsFactors = FALSE
  )
  
  membership_df <- data.frame(
    Gene.Symbol   = names(memb),
    Module = as.integer(memb),
    stringsAsFactors = FALSE
  )
  
  # Degree inside module only
  membership_df$Degree_inModule <- vapply(seq_along(memb), function(i) {
    node <- names(memb)[i]
    mod  <- memb[[i]]
    vs   <- names(memb[memb == mod])
    subg <- induced_subgraph(g, vids = vs)
    as.integer(degree(subg)[node])
  }, integer(1))
  
  return(list(summary_df = summary_df,
              membership_df = membership_df,
              edges_by_module = edges_by_module))
}


