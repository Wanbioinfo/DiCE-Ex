# ---------------------- To ensure consistent column naming for gene column -------------------------
# Ensure consistent column naming for gene column

#' Run normalize_dge_cols: Standardizes common differential-expression output columns to a fixed
#' schema used by DiCE.
#' Helper function (not for users)
#'
#' @param dge_data Dataframe of differential gene expression data
#'
#' @return Dataframe of differential gene expression data with renamed columns
#' @noRd
normalize_dge_cols <- function(dge_data){
  
  # Define possible matches
  gene_cols <- c("gene","genes","gene_name","gene_names","genename", "genenames","gene.name","gene.names",
                 "gene_symbol", "genesymbol", "gene.symbol", "gene_symbols", "genesymbols", "gene.symbols",
                 "symbol", "symbols")
  pval_cols <- c("pvalue","p.value","p_value",
                 "pval","p.val","p_val",
                 "pvals","p.vals","p_vals")
  adjP_cols <- c("adj.p.val","adj_pval","adjp","padj",
                 "padjusted","adjpval","adjpv",
                 "qval","q.value","qvalue","q_value",
                 "fdr","fdr_pval","p_val_adj")
  logfc_cols <- c("logfc","log2foldchange","logfc.value",
             "log2fc","log_fold_change","log2_fold_change",
             "log2ratio","log2_ratio","lfc",
             "avg_log2fc", "avg_logfc")
  
  colnames(dge_data)[tolower(colnames(dge_data)) %in% tolower(gene_cols)] <- "Gene.Symbol"
  colnames(dge_data)[tolower(colnames(dge_data)) %in% tolower(pval_cols)] <- "P.Value"
  colnames(dge_data)[tolower(colnames(dge_data)) %in% tolower(adjP_cols)] <- "adj.P.Val"
  colnames(dge_data)[tolower(colnames(dge_data)) %in% tolower(logfc_cols)] <- "logFC"
  
  return(dge_data)
}

