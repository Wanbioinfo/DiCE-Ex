#' Run DiCE Pipeline
#'
#' This function runs the DiCE pipeline for gene biomarker identification from transcriptomic data using differential expression, feature selection, and network-based centrality metrics.
#'
#' @param data_type Sequencing data type ("bulkRNA-seq" OR "scRNA-seq"). Default "bulkRNA-seq".
#' @param species "human" or "mouse". Default "human".
#' @param dge_file_path File path to the differential gene expression Rds/xlsx/csv/tsv file. The columns must include: "Gene.Symbol", "logFC", "P.Value", "adj.P.Val".
#' @param normGeneExp_file_path File path for the normalized gene expression Rds/xlsx/csv/tsv file (genes x samples/cells; make sure genes are under "Gene.Symbol" column).
#' @param rawGeneExp_file_path File path for the raw UMI counts Rds/xlsx/csv/tsv file (only needed when data_type = "scRNA-seq"). (genes x samples/cells; make sure genes are under "Gene.Symbol" column)
#' @param metadata_file_path Path to sample/cell metadata file (Rds/xlsx/csv/tsv). Must contain "Sample" or "Cell" and 'Phenotype' columns.
#' @param treatment Label of treatment samples (e.g., "Tumor").
#' @param control Label of control samples (e.g., "Normal").
#' @param remove_pc_genes If TRUE, protein-coding genes are removed from the input data before analysis. Default TRUE.
#' @param loose_criteria Statistical significance metric used for initial gene filtering ("P.Value" or "adj.P.Val"). Default "adj.P.Val".
#' @param loose_cutoff Numeric threshold for filtering based on 'loose_criteria'. Genes with values ≤ this cutoff are retained. Default 0.05.
#' @param logFC_cutoff Minimum absolute log2 fold change threshold for retaining genes. Default 0.
#' @param ig_method Character string specifying the Information Gain strategy to use. Options are:
#'   \describe{
#'     \item "IG" – Standard information gain
#'     \item "wIG" – Weighted information gain
#'     \item "none" – Skip information gain filtering
#'   }
#'   Default is "IG".
#' @param B Number of bootstrap resamples used in weighted IG calculation. Default is 300.
#' @param ig_cutoff Method for selecting IG-filtered genes. Options include:
#'   \describe{
#'     \item{"all_mean"}{Retain all genes with IG greater than the mean IG computed across all genes (including zeros). This is the standard mean-based filter. (Default)}
#'     \item{"all_median"}{Retain all genes with IG greater than the median IG computed across all genes (including zeros). More robust than the mean for extremely skewed IG distributions.}
#'     \item{"nonzero_mean"}{Compute the mean IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero mean threshold. Useful when IG has many exact zeros.}
#'     \item{"nonzero_median"}{Compute the median IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero median threshold. Ideal when IG is sparse and strongly right-skewed.}
#'     \item{"all_nonzero"}{Retain all genes with IG > 0. This keeps every gene carrying information with respect to the class label and excludes only those with zero information gain.}
#'     \item{"custom"}{Retain genes with IG greater than a user-defined threshold specified in \code{ig_custom_cutoff}}.
#'   }
#' @param ig_custom_cutoff Numeric IG threshold used only when \code{ig_cutoff = "custom"}.
#' @param corr_mode Mode for computing gene–gene correlation. Options are:
#' \describe{
#'   \item{"directCorr"}{Use raw or normalized expression values without dropping zero-expression cells. (Default)}
#'   \item{"remove_zerocells"}{Exclude cell pairs where both genes have zero expression before computing correlation.}
#'   \item{"zinbwave"}{Apply Zinb-Wave model for gene expression denoising before computing correlation (for scRNA-seq data).}
#'   \item{"newwave"}{Apply NewWave model for gene expression denoising before computing correlation (for scRNA-seq data).}
#' }
#' @param n_denoise_latent Number of latent factors used for single-cell normalization (K in zinbwave/newwave models). Default is 5.
#' @param max_denoise_iters Maximum number of iterations for the optimization step in denosing models (maxiter_optimize in zinbwave/newwave models). Default is 100).
#' @param corr_method Correlation method ("pearson" OR "spearman"). Default "pearson".
#' @param corr_pval_cutoff P-value threshold (default 1, no filtering); correlations with p-value > cutoff are set to 0.
#' @param ppi_db Database name for the PPI network. "stringdb" or "biogrid". Default is "stringdb".
#' @param stringDB_confidence Cutoff for confidence score in Stringdb PPI (Default is 400).
#' @param centrality_list Character vector of centrality metrics to compute.  
#'   Valid options include: "betweenness", "eigen vector", "pagerank", "closeness", "harmonic", "authority", "strength". 
#'   (Names are normalized internally.)
#' @param dice_rules A list of DiCE selection rules used to classify genes as final DiCE genes after ensemble ranking.
#'   
#'   Each rule must be created using helper functions such as
#'   \code{dice_centrality_rule()} or \code{dice_ensemble_rule()}.
#'   
#'   For a centrality rule, a gene passes if its centrality value is within the
#'   top specified cutoff in either the treatment network or the control network.
#'   For an ensemble rule, a gene passes if its \code{Ensemble_Rank} is within
#'   the specified cutoff.
#'   
#'   Example:
#'   \preformatted{
#'   dice_rules = list(
#'     dice_centrality_rule(
#'       metric = "betweenness",
#'       threshold_type = "percent",
#'       threshold = 25
#'     ),
#'     dice_ensemble_rule(
#'       threshold_type = "rank",
#'       threshold = 200
#'     )
#'   )
#'   }
#'
#' @param dice_logic Character string specifying how multiple \code{dice_rules}
#'   are combined to determine final DiCE genes. Must be either \code{"AND"} or
#'   \code{"OR"}.
#'   
#'   If \code{"AND"}, a gene must satisfy all rules.
#'   If \code{"OR"}, a gene must satisfy at least one rule.
#'   
#'   If only one rule is provided in \code{dice_rules}, \code{dice_logic} has no
#'   practical effect.
#' @param markTF If TRUE, transcription factor genes will be marked. Default FALSE.
#' 
#'
#' @return A named list containing:
#' \describe{
#'   \item{dice_results_df}{A data frame of final DiCE-ranked genes, including: Network centrality values across selected metrics, Ensemble ranking score,
#'   Information Gain (IG) or Weighted Information Gain (WIG), if applicable, Differential expression statistics (e.g., logFC, p-value, adjusted p-value),
#'   Transcription factor gene annotation, if applicable.}
#'   \item{interactions_df}{
#'   A data frame of protein–protein interactions used to construct the network,
#'   including interaction pairs and their corresponding edge weights
#'   (absolute correlation coefficients).
#'   }
#'   \item{phase2_denoised_geneExp}{
#'   Returned only when \code{data_type = "scRNA-seq"}.
#'   A filtered denoised gene expression matrix used for Phase II
#'   network construction and correlation analysis.
#'   }
#' }
#' 
#' @importFrom dplyr mutate filter as_tibble bind_rows select transmute group_by ungroup rename
#' @importFrom stringr str_to_lower str_detect
#' @importFrom FSelectorRcpp information_gain
#' @importFrom igraph graph_from_data_frame set_edge_attr E E<- V V<-
#' @importFrom igraph eigen_centrality betweenness authority_score closeness
#' @importFrom igraph harmonic_centrality page_rank strength simplify
#' @importFrom igraph cluster_louvain membership sizes modularity cluster_leiden
#' @importFrom igraph induced_subgraph degree vcount ecount as_adj
#' @importFrom utils stack combn read.csv read.delim read.table
#' @importFrom parallel mclapply
#' @importFrom NetWeaver ensemble_rank
#' @importFrom scuttle logNormCounts
#' @importFrom stats cor quantile median as.formula setNames sd
#' @importFrom BiocParallel register MulticoreParam
#' @importFrom readxl read_excel
#' @importFrom data.table fread
#' @importFrom AnnotationDbi keys
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom tidyr unnest_wider
#' @importFrom Matrix rowSums sparseMatrix
#' @importFrom Hmisc rcorr
#'
#' @examples
#' \dontrun{
#' dice_results <- perform_DiCE(
#'   data_type = "bulkRNA-seq",
#'   dge_file_path = "path/to/dge_results.rds",
#'   normGeneExp_file_path = "path/to/logNorm_geneExp.rds",
#'   metadata_file_path = "path/to/metadata.rds",
#'   treatment = "Tumor",
#'   control = "Normal",
#'   loose_criteria = "adj.P.Val",
#'   loose_cutoff = 0.05,
#'   logFC_cutoff = 1,
#'   species = "human",
#'   ig_method = "IG",
#'   B = 200,
#'   ig_cutoff = "nonzero_mean",
#'   corr_mode = "directCorr",
#'   corr_method = "pearson",
#'   corr_pval_cutoff = 0.05,
#' 
#' )
#' head(dice_results$dice_results_df)
#' }
#'
#' @export
perform_DiCE <- function(
    data_type = "bulkRNA-seq",
    species = "human",
    dge_file_path = NULL,
    normGeneExp_file_path = NULL,
    rawGeneExp_file_path = NULL,
    metadata_file_path = NULL,
    treatment = NULL,
    control = NULL,
    remove_pc_genes = TRUE,
    loose_criteria = "adj.P.Val",
    loose_cutoff = 0.05,
    logFC_cutoff = 0,
    ig_method = "IG",
    B = 300,
    ig_cutoff = "all_mean",
    ig_custom_cutoff = NULL, 
    corr_mode = "directCorr",
    n_denoise_latent = 5,
    max_denoise_iters = 100,
    corr_method = "pearson",
    corr_pval_cutoff = 1,
    ppi_db = "stringdb",
    stringDB_confidence = 400,
    centrality_list = c("betweenness", "eigen vector"),
    dice_rules = list(
      dice_centrality_rule(
        metric = "betweenness",
        threshold_type = "percent",
        threshold = 25
      ),
      dice_centrality_rule(
        metric = "eigen vector",
        threshold_type = "percent",
        threshold = 25
      )
    ),
    dice_logic = "AND",
    markTF = FALSE
    ) 
{
  
  dge_data <- NULL
  normGeneExp_data <- NULL
  rawGeneExp_data <- NULL


  # String db downloaded files
  if(tolower(species) == "human"){
    string_protInfo_file <- "extdata/stringDB_v12/human/9606.protein.info.v12.0.txt"
    string_ppi_file <- "extdata/stringDB_v12/human/9606.protein.links.v12.0.txt.gz"
    
    biogrid_ppi_file <- "extdata/biogrid/human/BIOGRID-ORGANISM-Homo_sapiens-5.0.257.tab.txt.gz"
    
    taxonID <- 9606
    
  }else if(tolower(species) == "mouse"){
    string_protInfo_file <- "extdata/stringDB_v12/mouse/10090.protein.info.v12.0.txt"
    string_ppi_file <- "extdata/stringDB_v12/mouse/10090.protein.links.v12.0.txt.gz"
    
    biogrid_ppi_file <- "extdata/biogrid/mouse/BIOGRID-ORGANISM-Mus_musculus-5.0.257.tab.txt.gz"
    
    taxonID = 10090
    
  }else{
    stop("Invalid species!. DiCE supports only for 'human' and 'mouse'")
  }

  ############################### Load data #############################################
  
  message("Loading input data...")
  
  # Load DGE data and gene expression data
  if (is.null(dge_file_path) || !file.exists(dge_file_path)) {
    stop("Invalid or missing 'dge_file_path'")
  }
  
  dge_data <- read_any(dge_file_path)
  
  
  # Load meta data file
  if (is.null(metadata_file_path) || !file.exists(metadata_file_path)) {
    stop("Invalid or missing 'metadata_file_path'")
  }
  
  meta_data <- read_any(metadata_file_path)
  
  if (data_type == "scRNA-seq") {
    
    
    # read raw read counts
    if (is.null(rawGeneExp_file_path) || !file.exists(rawGeneExp_file_path)) {
      stop("Invalid or missing 'rawGeneExp_file_path' for scRNA-seq data.")
    }
    
    rawGeneExp_data <- read_any(rawGeneExp_file_path)
    
    
    # read normalized 
    if (tolower(ig_method) %in% c("ig", "wig")) {
      
      if (is.null(normGeneExp_file_path) || !file.exists(normGeneExp_file_path)) {
        stop("Invalid or missing 'normGeneExp_file_path' when ig_method is 'ig' or 'wig'.")
      }
      
      normGeneExp_data <- read_any(normGeneExp_file_path)
      
    }
    
    if (!(corr_mode %in% c("ZINB-WaVE", "NewWave"))) {
      
      if (is.null(normGeneExp_file_path) || !file.exists(normGeneExp_file_path)) {
        stop("Invalid or missing 'normGeneExp_file_path' when ig_method is 'ig' or 'wig'.")
      }
      
      normGeneExp_data <- read_any(normGeneExp_file_path)

    }
    
  } else if (data_type == "bulkRNA-seq") {
    if (is.null(normGeneExp_file_path) || !file.exists(normGeneExp_file_path)) {
      stop("Invalid or missing 'normGeneExp_file_path'.")
    }
    normGeneExp_data <- read_any(normGeneExp_file_path)
  
    
  } else {
    stop("Invalid value for 'data_type'. Select from 'bulkRNA-seq' or 'scRNA-seq'.")
  }
  
  # Change column names to DiCE col names
  dge_data <- normalize_dge_cols(dge_data)
  
  ###################### Remake the gene expression data ################################
  
  if (!is.null(rawGeneExp_data)){
    rawGeneExp_data <- exp_to_DiCE_format(meta_data,rawGeneExp_data)
    rawGeneExp_data <- NA_to_numeric(rawGeneExp_data)
  }
  
  if (!is.null(normGeneExp_data)){
    normGeneExp_data <- exp_to_DiCE_format(meta_data,normGeneExp_data)
    normGeneExp_data <- NA_to_numeric(normGeneExp_data)
  }
  
  
  
  
  
  
  ###################### Keep only protein coding genes #################################
  
  if (isTRUE(remove_pc_genes) || tolower(remove_pc_genes) %in% c("yes", "true", "y", "1")){
    message("Starting protein-coding gene filtering...")
    
    dge_data <- keep_protCoding_dgeData(taxonID, dge_data)
    
    prot_genes <- dge_data$Gene.Symbol
    
    if (!is.null(normGeneExp_data) && nrow(normGeneExp_data) > 0) {
      normGeneExp_data <- keep_protCoding_expData(prot_genes, normGeneExp_data)
    }
    
    if (data_type == "scRNA-seq") {
      rawGeneExp_data <- keep_protCoding_expData(prot_genes, rawGeneExp_data)
    }
    
  }
  
  ############################### Phase 1 ###############################################
  
  # Construction of a candidate gene pool by Differential Gene Expression
  # analysis with a loose cutoff
  
  phase1_res <- run_phase1(dge_data, loose_criteria, loose_cutoff, logFC_cutoff)
  
  phase1_criteria <- paste0(loose_criteria,
                            "<",
                            loose_cutoff,
                            " & ",
                            "|logFC|>",
                            logFC_cutoff)
  
  message(paste0("#Genes in Phase1 (", phase1_criteria, ") = ", nrow(phase1_res)))
  
  
  ############################### Phase 2 ###############################################
  
  # Selection of the top discriminative genes using the Information Gain
  # filter approach and construct the PPI interactions.
  if (tolower(ig_method) %in% c("ig", "wig")){
    phase2_res <- run_phase2(normGeneExp_data, phase1_res, ig_method, ig_cutoff, ig_custom_cutoff, B)
    
    phase2_genes_df <- phase2_res$phase2_genes_df
    ig_df <- phase2_res$infoGain_df
    
    phase1_res <- merge(phase1_res, ig_df, by = "Gene.Symbol", all.x = TRUE) 
    
    if (tolower(ig_method) == "wig") {
      ig_term <- "WeightedIG"
    }else{
      ig_term <- "IG"
    }
    
    if (ig_cutoff == "all_mean"){
      phase2_criteria <- paste0("Mean(",ig_term,")")
      
    }else if (ig_cutoff == "all_median"){
      phase2_criteria <- paste0("Median(",ig_term,")")
      
    }else if (ig_cutoff == "nonzero_mean"){
      phase2_criteria <- paste0("Mean(",ig_term,">0)")
      
    }else if (ig_cutoff == "nonzero_median"){
      phase2_criteria <- paste0("Median(",ig_term,">0)")
    
    }else if (ig_cutoff == "all_nonzero"){
      phase2_criteria <- paste0(ig_term,">0")
      
    }else if (ig_cutoff == "custom"){
      phase2_criteria <- paste0(ig_term,">",ig_custom_cutoff)
      
    }else{
      stop("Invalid cutoff for IG!. Please select from 'all_mean', 'all_median', 'nonzero_mean', 'nonzero_median', 'all_nonzero', and 'custom'.")
    }
    
    message(paste0("#Genes in Phase2 (",phase2_criteria,") = " , nrow(phase2_genes_df)))
    
    to_phase3 <- phase2_genes_df
    
  }else if (tolower(ig_method) == "none") {
    message("#Phase 2 skipped: ig_method = 'none'")
    phase2_genes_df <- data.frame()
    to_phase3 <- phase1_res
    
  }else{
    stop("Invalid value for 'ig_method'. Allowed values are: 'IG', 'wIG', or 'none'.")
  }
  
  ############################### Phase 3 ###############################################
  
  # Making weighted PPI using (1-|C.C|) for each phenotype
  
  
  
  # Extract interactions
  if(tolower(ppi_db) == "stringdb"){
    message("Phase3: Creating PPI network for Phase2 genes from StringDB")
    
    ppi_results_phase3 <- extract_stringdb_PPI(string_protInfo_file, 
                                               string_ppi_file, 
                                               stringDB_confidence,
                                               to_phase3)
  }else if (tolower(ppi_db) == "biogrid"){
    message("Phase3: Creating PPI network for Phase2 genes from BioGRID")
    
    ppi_results_phase3 <- extract_biogrid_PPI(biogrid_ppi_file,
                                              taxonID,
                                              to_phase3)
    
  }else{
    stop("Invalid value for 'ppi_db'. Allowed values are: 'stringdb', or 'biogrid'.")
  }
  

  pp_interactions <- ppi_results_phase3$interactions
  mapped_proteins <- ppi_results_phase3$mapped_proteins
  
  
  if (!is.null(normGeneExp_data) && nrow(normGeneExp_data) > 0) {
    filtered_normGeneExp_data <- getNetworkGene_expression(normGeneExp_data,
                                                           mapped_proteins)
  }
    
  

  if (data_type == "scRNA-seq"){
    filtered_raw_geneExp <- getNetworkGene_expression(rawGeneExp_data,
                                                      mapped_proteins)
  }
  
  if (corr_mode == "zinbwave"){
    print("Raw gene expression - zinbwave denoising")
    filtered_denoised_geneExp <- zinbWave_model(raw_geneExp = filtered_raw_geneExp, 
                                                K = n_denoise_latent,
                                                maxIters = max_denoise_iters)
    geneExp_to_corr <- filtered_denoised_geneExp
    
  }else if (corr_mode == "newwave"){
    print("Raw gene expression - NewWaVE denoising")
    filtered_denoised_geneExp <- newWave_model(raw_geneExp = filtered_raw_geneExp, 
                                               K = n_denoise_latent,
                                               maxIters = max_denoise_iters)
    geneExp_to_corr <- filtered_denoised_geneExp
    
  }else if ((corr_mode == "directCorr") | (corr_mode == "remove_zerocells")){
    print("Normalized gene expressions")
    geneExp_to_corr <- filtered_normGeneExp_data
    
  }else{
    stop("Invalid correlation mode. Use 'directCorr', 'remove_zerocells', 'zinbwave', or 'newwave'")
  }
  
  # Run Phase 3
  if (is.null(treatment)){
    stop("Missing 'treatment'")
  }
  
  if (is.null(control)){
    stop("Missing 'control'")
  }
  
  
  phase3_interactions_df <- run_phase3(pp_interactions, geneExp_to_corr, treatment, control, 
                           corr_method, corr_mode, corr_pval_cutoff)
  
  phase3_genes <- unique(c(phase3_interactions_df$source, phase3_interactions_df$target))

  if(tolower(ppi_db) == "stringdb"){
    message(paste0("#Genes in Phase3 (StringDB PPI - ",stringDB_confidence," ) = "),length(phase3_genes))
  }else{
    message(paste0("#Genes in Phase3 (BioGRID PPI) = "),length(phase3_genes))
  }

  
  ############################### Phase 4 ###############################################
  
  # Calculate network centralities
  
  phase4_centralities_df <- run_phase4(phase3_interactions_df, centrality_list)
  
  phase4_centralities_df$Phase <- "III/IV"
  
  phase4_centralities_out_df <- merge(phase4_centralities_df, mapped_proteins, 
                                      by.x = "Gene.Symbol", by.y = "Gene.Symbol")
  colnames(phase4_centralities_out_df)[1] <- "Gene.Symbol"
  
  # Rearrange the columns
  base_cols <- c("Gene.Symbol", 
                 if ("STRING_id" %in% colnames(phase4_centralities_out_df)) "STRING_id",
                 if ("BioGRID_id" %in% colnames(phase4_centralities_out_df)) "BioGRID_id",
                 "logFC", 
                 "adj.P.Val", 
                 "P.Value", 
                 if ("IG" %in% colnames(phase4_centralities_out_df)) "IG")
  cent_cols <- colnames(phase4_centralities_out_df)[2:((length(centrality_list)*4)+1)]
  tail_cols <- "Phase"
  final_cols <- c(base_cols, cent_cols, tail_cols)
  
  phase4_centralities_out_df <- phase4_centralities_out_df[, final_cols, drop = FALSE]
  
  message(paste0("#Genes in Phase4 (Centralities/Ensemble ranking) = ", nrow(phase4_centralities_df)))
  
  # Ensemble ranking
  dice_results_df <- ensemble_Ranking(dge_data, phase1_res, 
                                      phase2_genes_df, 
                                      phase4_centralities_out_df,
                                      centrality_list)
  

  ############################### Phase 5 ###############################################
  
  dice_results_df <- apply_dice_rules(
    dice_results_df = dice_results_df,
    dice_rules = dice_rules,
    dice_logic = dice_logic
  )
  
  dice_genes <- dice_results_df[dice_results_df$Phase == "DiCE",]$Gene.Symbol
  rule_text <- format_dice_rules_text(dice_rules, dice_logic)
  message(paste0("#Genes in Phase5 (", rule_text, ") = ", length(dice_genes)))

  # rename IG if the ig_method is weighted information gain
  if (tolower(ig_method) == "wig" && "IG" %in% colnames(dice_results_df)) {
    colnames(dice_results_df)[colnames(dice_results_df) == "IG"] <- "WIG"
  }
  
  # rename treatment and control column names to original names
  colnames(dice_results_df) <- gsub("treatment", treatment, colnames(dice_results_df))
  colnames(dice_results_df) <- gsub("control", control, colnames(dice_results_df))
  
  ######################## Annotate Transcription factors ###############################
  
  
  if (isTRUE(markTF) || tolower(markTF) %in% c("yes", "true", "y", "1")){
    dice_results_df <- annotateTFs(dice_results_df, tf_file)
  }
  
  
  ############################### Return output ########################################
  
  # remove PPI id col
  # Remove STRING_id and BioGRID_id columns if they exist
  cols_to_remove <- c("STRING_id", "BioGRID_id")
  dice_results_df <- dice_results_df[, !colnames(dice_results_df) %in% cols_to_remove]
  
  # rename treatment and control column names to original names
  colnames(phase3_interactions_df) <- gsub("treatment", treatment, colnames(phase3_interactions_df))
  colnames(phase3_interactions_df) <- gsub("control", control, colnames(phase3_interactions_df))
  
  if (data_type == "scRNA-seq" & tolower(corr_mode) %in% c("zinbwave", "newwave")) {
    return(list(
      dice_results_df = dice_results_df,
      interactions_df = phase3_interactions_df,
      phase2_denoised_geneExp = filtered_denoised_geneExp
    ))
  }
  
  return(list(dice_results_df = dice_results_df,
         interactions_df = phase3_interactions_df))

}
