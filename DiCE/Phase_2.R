# Phase II: Selection of the top discriminative genes using the Information Gain
# filter approach.


#' Run Phase 2: Gene filtering based on information gain
#' Helper function (not for users)
#'
#' @param geneExp_data Dataframe of gene expression data for each sample with the sample group.
#' @param phase1_res DiCE phase1 results.
#' @param ig_method Character string specifying the Information Gain strategy to use. IG or wIG or none. Default is IG.
#' @param ig_cutoff Method for selecting IG-filtered genes.
#' @param wrs_B Number of bootstrap resamples used in weighted IG calculation. Default is 300.
#' @param wrs_seed Integer seed used to control randomness.
#'
#' @return Filtered dataframe of Phase2 candidate genes with their IG values.
#' @noRd
run_phase2 <- function(geneExp_data, phase1_res, ig_method, ig_cutoff, ig_custom_cutoff = NULL,
                       wrs_B = 300, wrs_seed = 123){
  
  # inputs and basic checks
  if(ncol(geneExp_data) < 2){
    stop("Gene expression data must have genes + a class column (last column).")
  }
  
  class <- geneExp_data[[ncol(geneExp_data)]]
  
  if(is.numeric(class)){
    stop("Please check your last column in the gene expression data. ",
         "It should represent group or class of the samples and be categorical (factor).",
         call. = FALSE)
  }
  
  class <- as.factor(class)
  
  if (length(levels(class))!=2){
    stop("Please check your last column in the gene expression data. ",
         "It should contain only two group or class levels.",
         call. = FALSE)
  }
  
  # Filter the genes expression from Phase1
  phase1_genes <- phase1_res$Gene.Symbol
  
  if (length(phase1_genes) == 0) stop("No Phase I genes found in geneExp_data columns.")
  
  geneExp_data <- geneExp_data[,colnames(geneExp_data) %in% phase1_genes]
  rownames(geneExp_data) <- NULL
  new_geneExp_data <- cbind(as.data.frame(geneExp_data), class)
  
  # Calculate the information gain with warning handler
  # weighted or traditional IG
  if (tolower(ig_method) == "wig") {
    # Weighted resampling IG (approximates class-weighted IG; keeps MDL discretization)
    wrs_tbl <- IG_weighted_resample_fast_repro(
      df = new_geneExp_data,
      class_col = "class",
      B = wrs_B,
      seed = wrs_seed,
      n_cores_boot = 8   # use ~all cores
    )
    
    # Use median IG as the score (robust), rename to IG for downstream compatibility
    infoGain_ordered_df <- wrs_tbl %>%
      transmute(Gene.Symbol,
                IG = IG_wrs_median) %>%
      arrange(desc(IG))
    
  } else {
    # Traditional IG (MDL discretization inside FSelectorRcpp)
    infoGain_df <- information_gain(class ~ ., new_geneExp_data, type = "infogain")
    
    rownames(infoGain_df) <- infoGain_df$attributes
    infoGain_df <- infoGain_df[, "importance", drop = FALSE]
    colnames(infoGain_df)[1] <- "IG"
    
    # Order by IG (descending)
    infoGain_ordered_df <- infoGain_df %>%
      tibble::rownames_to_column("Gene.Symbol") %>%
      arrange(desc(IG))
  }
    
  
  # Ranks and Phase II selection (>= mean IG)
  gene_rankings <- infoGain_ordered_df %>%
    mutate(Rank = seq_len(n()))
  
  ig_cutoff <- tolower(ig_cutoff)
  
  # ig cutoff all_mean, all_median, nonzero_mean, nonzero_median
  if (ig_cutoff == "all_mean") {
    
    cutoff_genes <- gene_rankings %>%
      filter(IG >= mean(IG, na.rm = TRUE))
    
  } else if (ig_cutoff == "all_median") {
    
    cutoff_genes <- gene_rankings %>%
      filter(IG >= median(IG, na.rm = TRUE))
    
  } else if (ig_cutoff == "nonzero_mean") {
    
    cutoff_genes <- gene_rankings %>%
      filter(IG >= mean(IG[IG > 0], na.rm = TRUE))
    
  } else if (ig_cutoff == "nonzero_median") {
    
    cutoff_genes <- gene_rankings %>%
      filter(IG >= median(IG[IG > 0], na.rm = TRUE))
    
  } else if (ig_cutoff == "all_nonzero") {
    
    cutoff_genes <- gene_rankings %>%
      filter(IG > 0)
    
  } else if (ig_cutoff == "custom") {
    
    if (is.null(ig_custom_cutoff) || 
        !is.numeric(ig_custom_cutoff) || 
        length(ig_custom_cutoff) != 1) {
      stop("When ig_cutoff = 'custom', ig_custom_cutoff must be a single numeric value.")
    }
    
    cutoff_genes <- gene_rankings %>%
      filter(IG > ig_custom_cutoff)
    
  } else {
    
    stop("Invalid cutoff for IG!. Please select from 'all_mean', 'all_median', 'nonzero_mean', 'nonzero_median', 'all_nonzero', and 'custom'.")
  }
  
  
  phase2_genes_df <- merge(cutoff_genes, phase1_res, by = "Gene.Symbol")
  phase2_genes_df <- phase2_genes_df[order(phase2_genes_df$Rank), ]
  phase2_genes_df$Phase <- "II"
  
  return(list(
    phase2_genes_df = phase2_genes_df,
    infoGain_df     = infoGain_ordered_df   # columns: Gene.Symbol, IG
  ))
  
}
