#!/usr/bin/env Rscript

message("===================================")
message("DiCE CLI Run")
message("Start time: ", Sys.time())
message("===================================")

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(openxlsx)
  library(dplyr)
  library("tibble")
  library("stringr")
  library('tidyr')
  library('purrr')
  library('Matrix')
  library("FSelectorRcpp")
  library("igraph")
  library("data.table")
  library("NetWeaver")
  library("praznik")
  library("reticulate")
  library("stats")
  library("utils")
  library("parallel")
  library("annotate")
  library("org.Hs.eg.db")
  library("AnnotationDbi")
  library("org.Mm.eg.db")
  library("readxl")
  library('Hmisc')
})

# ---------- Helpers ----------
die <- function(...) {
  msg <- paste0(...)
  message("ERROR: ", msg)
  quit(status = 1)
}

has_value <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && nzchar(x)

norm_path <- function(x) {
  if (!has_value(x)) return(NULL)
  normalizePath(path.expand(x), mustWork = FALSE)
}

has_file <- function(x) has_value(x) && file.exists(x)

ensure_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

parse_dice_rules <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1) {
    if (is.na(x) || !nzchar(trimws(x))) return(NULL)
    x <- tryCatch(
      jsonlite::fromJSON(x, simplifyVector = FALSE),
      error = function(e) die("Failed to parse --dice_rules JSON: ", e$message)
    )
  }
  if (is.data.frame(x)) {
    x <- split(x, seq_len(nrow(x)))
    x <- lapply(x, function(row) as.list(row))
  }
  if (!is.list(x) || length(x) == 0) {
    die("dice_rules must be a non-empty list or a valid JSON string.")
  }
  x
}

normalize_species <- function(x) {
  x <- tolower(trimws(x))
  if (x %in% c("human", "hs", "homo_sapiens", "homo sapiens")) return("human")
  if (x %in% c("mouse", "mm", "mus_musculus", "mus musculus")) return("mouse")
  die("Unsupported species: ", x, " (use 'human' or 'mouse').")
}

# Helper to safely extract nested module config with defaults
get_mod_cfg <- function(modules_cfg, branch, key, default) {
  val <- modules_cfg[[branch]][[key]]
  if (is.null(val)) default else val
}

# ---------- CLI options ----------
option_list <- list(
  make_option(c("--species"),                type="character", default="human",          dest="species"),
  make_option(c("--dge_file_path"),          type="character", default=NULL,             dest="dge_file_path"),
  make_option(c("--normGeneExp_file_path"),  type="character", default=NULL,             dest="normGeneExp_file_path"),
  make_option(c("--metadata_file_path"),     type="character", default=NULL,             dest="metadata_file_path"),
  make_option(c("--treatment"),              type="character", default=NULL,             dest="treatment"),
  make_option(c("--control"),                type="character", default=NULL,             dest="control"),
  make_option(c("--loose_criteria"),         type="character", default="adj.P.Val",      dest="loose_criteria"),
  make_option(c("--loose_cutoff"),           type="double",    default=0.05,             dest="loose_cutoff"),
  make_option(c("--logFC_cutoff"),           type="double",    default=0,                dest="logFC_cutoff"),
  make_option(c("--ig_method"),              type="character", default="IG",             dest="ig_method"),
  make_option(c("--B"),                      type="integer",   default=300,              dest="B"),
  make_option(c("--ig_cutoff"),              type="character", default="all_mean",        dest="ig_cutoff"),
  make_option(c("--ig_custom_cutoff"),       type="double",    default=NA,               dest="ig_custom_cutoff"),
  make_option(c("--corr_method"),            type="character", default="pearson",         dest="corr_method"),
  make_option(c("--centrality_list"),        type="character", default="betweenness,eigenvector", dest="centrality_list"),
  make_option(c("--dice_rules"),             type="character", default=NULL,             dest="dice_rules"),
  make_option(c("--dice_logic"),             type="character", default="AND",            dest="dice_logic"),
  make_option(c("--seed"),                   type="integer",   default=123,              dest="seed"),
  make_option(c("--outdir"),                 type="character", default="dice_cli_output", dest="outdir"),
  make_option(c("--job_name"),               type="character", default="DiCE_CLI_Run",   dest="job_name"),
  make_option(c("--config"),                 type="character", default=NULL,             dest="config"),
  make_option(c("--ppi_db"),                 type="character", default="stringdb",       dest="ppi_db",
              help="PPI database to use: stringdb or biogrid [default %default]"),
  make_option(c("--remove_pc_genes"),        type="logical",   default=TRUE,            dest="remove_pc_genes",
              help="Remove protein-coding genes filter [default %default]"),
  make_option(c("--stringDB_confidence"),    type="integer",   default=400,             dest="stringDB_confidence",
              help="StringDB confidence score threshold (0-1000) [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------- Load config JSON ----------
if (!is.null(opt$config)) {
  if (!file.exists(opt$config)) die("Config not found: ", opt$config)
  cfg <- jsonlite::fromJSON(opt$config, simplifyVector = FALSE)
  message("Config keys: ", paste(names(cfg), collapse = ", "))
  for (nm in names(cfg)) opt[[nm]] <- cfg[[nm]]
}

# ---------- Resolve paths ----------
opt$normGeneExp_file_path <- norm_path(opt$normGeneExp_file_path)
opt$dge_file_path         <- norm_path(opt$dge_file_path)
opt$metadata_file_path    <- norm_path(opt$metadata_file_path) 

# ---------- Validate required fields ----------
if (!has_file(opt$normGeneExp_file_path)) die("--normGeneExp_file_path missing or file not found.")
if (!has_file(opt$dge_file_path))         die("--dge_file_path missing or file not found.")
if (!has_file(opt$metadata_file_path)) { die("--metadata_file_path missing or file not found.")}

if (!has_value(opt$treatment))            die("--treatment is required.")
if (!has_value(opt$control))              die("--control is required.")

species <- normalize_species(opt$species)
set.seed(opt$seed)
ensure_dir(opt$outdir)

# ---------- Parse modules config ----------
# Supports new nested format:
# "modules": { "run": true, "unweighted": {...}, "weighted": {...} }
# Falls back to flat "run_modules": true for backward compatibility

modules_cfg <- opt$modules  # will be NULL if not in config

run_modules     <- FALSE
run_unweighted  <- FALSE
run_weighted    <- FALSE

unweighted_algo        <- "louvain"
unweighted_res         <- 1
unweighted_leiden_itrs <- 10
unweighted_leiden_beta <- 0.01

weighted_algo        <- "louvain"
weighted_res         <- 1
weighted_leiden_itrs <- 10
weighted_leiden_beta <- 0.01

if (!is.null(modules_cfg) && is.list(modules_cfg)) {
  # New nested format
  run_modules    <- isTRUE(modules_cfg$run)
  run_unweighted <- run_modules && isTRUE(get_mod_cfg(modules_cfg, "unweighted", "enabled", TRUE))
  run_weighted   <- run_modules && isTRUE(get_mod_cfg(modules_cfg, "weighted",   "enabled", TRUE))
  
  unweighted_algo        <- get_mod_cfg(modules_cfg, "unweighted", "algorithm",   "louvain")
  unweighted_res         <- get_mod_cfg(modules_cfg, "unweighted", "resolution",  1)
  unweighted_leiden_itrs <- get_mod_cfg(modules_cfg, "unweighted", "leiden_itrs", 10)
  unweighted_leiden_beta <- get_mod_cfg(modules_cfg, "unweighted", "leiden_beta", 0.01)
  
  weighted_algo        <- get_mod_cfg(modules_cfg, "weighted", "algorithm",   "louvain")
  weighted_res         <- get_mod_cfg(modules_cfg, "weighted", "resolution",  1)
  weighted_leiden_itrs <- get_mod_cfg(modules_cfg, "weighted", "leiden_itrs", 10)
  weighted_leiden_beta <- get_mod_cfg(modules_cfg, "weighted", "leiden_beta", 0.01)
  
} else if (isTRUE(opt$run_modules)) {
  # Backward-compatible flat format
  run_modules    <- TRUE
  run_unweighted <- TRUE
  run_weighted   <- TRUE
}

# ---------- Source DiCE ----------
repo_root <- getwd()
dice_dir  <- file.path(repo_root, "DiCE")
if (!dir.exists(dice_dir)) dice_dir <- "DiCE"
if (!dir.exists(dice_dir)) die("Could not find DiCE directory. Expected: ", dice_dir)

source(file.path(dice_dir, "arrange_Input.R"))
source(file.path(dice_dir, "calculate_NetCentralities.R"))
source(file.path(dice_dir, "calculate_weightedIG.R"))
source(file.path(dice_dir, "corr_calculations.R"))
source(file.path(dice_dir, "ensemble_Ranking.R"))
source(file.path(dice_dir, "createPPI.R"))
source(file.path(dice_dir, "DiCE.R"))
source(file.path(dice_dir, "fileReader.R"))
source(file.path(dice_dir, "normalize_df_cols.R"))
source(file.path(dice_dir, "Phase_1.R"))
source(file.path(dice_dir, "Phase_2.R"))
source(file.path(dice_dir, "Phase_3.R"))
source(file.path(dice_dir, "Phase_4.R"))
source(file.path(dice_dir, "Phase_5.R"))
source(file.path(dice_dir, "protein_coding_filter.R"))
source(file.path(dice_dir, "unweighted_module_analysis.R"))
source(file.path(dice_dir, "weighted_module_analysis_helpers.R"))
source(file.path(dice_dir, "weighted_module_analysis.R"))

# ---------- Parse centralities ----------
if (is.list(opt$centrality_list)) {
  centrality_list <- unlist(opt$centrality_list)
} else {
  centrality_list <- trimws(strsplit(opt$centrality_list, ",", fixed = TRUE)[[1]])
}
centrality_list <- centrality_list[nzchar(centrality_list)]
if (length(centrality_list) == 0) die("No centralities provided in --centrality_list")

# ---------- Parse dice_rules ----------
opt$dice_logic <- toupper(trimws(opt$dice_logic))
if (!(opt$dice_logic %in% c("AND", "OR"))) die("--dice_logic must be 'AND' or 'OR'.")

ig_custom_cutoff <- if (tolower(opt$ig_cutoff) == "custom") opt$ig_custom_cutoff else NULL
if (tolower(opt$ig_cutoff) == "custom" && (is.null(ig_custom_cutoff) || is.na(ig_custom_cutoff))) {
  die("ig_cutoff=custom but --ig_custom_cutoff is missing.")
}

dice_rules <- parse_dice_rules(opt$dice_rules)

# ---------- Log settings ----------
message("normGeneExp_file_path : ", opt$normGeneExp_file_path)
message("dge_file_path         : ", opt$dge_file_path)
message("metadata_file_path    : ", opt$metadata_file_path)
message("treat / control       : ", opt$treatment, " / ", opt$control)
message("species               : ", species)
message("centralities          : ", paste(centrality_list, collapse = ", "))
message("dice_logic            : ", opt$dice_logic)
message("run_unweighted modules: ", run_unweighted)
message("run_weighted modules  : ", run_weighted)
message("ppi_db                : ", opt$ppi_db)
message("stringDB_confidence   : ", opt$stringDB_confidence)
message("remove_pc_genes       : ", opt$remove_pc_genes)

# ---------- Run DiCE ----------
start_time <- proc.time()

dice_res <- perform_DiCE(
  data_type             = "bulkRNA-seq",
  species               = species,
  dge_file_path         = opt$dge_file_path,
  normGeneExp_file_path = opt$normGeneExp_file_path,
  metadata_file_path    = opt$metadata_file_path,
  treatment             = opt$treatment,
  control               = opt$control,
  remove_pc_genes       = isTRUE(opt$remove_pc_genes),
  loose_criteria        = opt$loose_criteria,
  loose_cutoff          = opt$loose_cutoff,
  logFC_cutoff          = opt$logFC_cutoff,
  ig_method             = opt$ig_method,
  B                     = opt$B,
  ig_cutoff             = opt$ig_cutoff,
  ig_custom_cutoff      = ig_custom_cutoff,
  corr_mode             = "directCorr",
  corr_method           = opt$corr_method,
  corr_pval_cutoff      = 1,
  ppi_db                = opt$ppi_db,
  stringDB_confidence   = as.integer(opt$stringDB_confidence),
  centrality_list       = centrality_list,
  dice_rules            = dice_rules,
  dice_logic            = opt$dice_logic
)

dice_runtime <- proc.time() - start_time
message("DiCE runtime: ", round(dice_runtime[["elapsed"]] / 60, 2), " min")

result_df           <- as.data.frame(dice_res$dice_results_df)
phase3_interactions <- as.data.frame(dice_res$interactions_df)

# ---------- Save DiCE results ----------
out_xlsx <- file.path(opt$outdir, paste0(opt$job_name, "_dice_results.xlsx"))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "DiCE_results")
openxlsx::writeData(wb, "DiCE_results", result_df)
openxlsx::addWorksheet(wb, "Phase3_interactions")
openxlsx::writeData(wb, "Phase3_interactions", phase3_interactions)
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved DiCE results: ", out_xlsx)

# ---------- Resolve DiCE genes for module detection ----------
if (run_unweighted || run_weighted) {
  
  phase_col <- if ("Phase" %in% names(result_df)) "Phase" else
    if ("phase" %in% names(result_df)) "phase" else NULL
  if (is.null(phase_col)) die("Could not find Phase column in DiCE results.")
  
  dice_genes_df <- result_df[result_df[[phase_col]] == "DiCE", , drop = FALSE]
  if (nrow(dice_genes_df) == 0) die("No rows with Phase == 'DiCE'; cannot run module detection.")
  
  gene_col <- c("Gene.Symbol", "Gene.Name", "Gene", "gene")
  gene_col <- gene_col[gene_col %in% names(dice_genes_df)][1]
  if (is.na(gene_col)) die("Could not find a gene symbol column in DiCE results.")
  
  dice_genes <- as.character(dice_genes_df[[gene_col]])
  message("DiCE genes for module detection: ", length(dice_genes))
}

# ---------- Unweighted modules ----------
if (run_unweighted) {
  message("Running unweighted module detection...")
  uw_start <- proc.time()
  
  modules <- detect_PPI_unweightedModules(
    gene_list    = dice_genes,
    interactions = phase3_interactions,
    algorithm    = unweighted_algo,
    resolution   = unweighted_res,
    leiden_itrs  = unweighted_leiden_itrs,
    leiden_beta  = unweighted_leiden_beta,
    seed         = opt$seed
  )
  
  uw_runtime <- proc.time() - uw_start
  message("Unweighted modules runtime: ", round(uw_runtime[["elapsed"]] / 60, 2), " min")
  
  mod_xlsx <- file.path(opt$outdir, paste0(opt$job_name, "_unweighted_modules.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Module summary")
  openxlsx::writeData(wb, "Module summary",       as.data.frame(modules$summary_df))
  openxlsx::addWorksheet(wb, "Module statistics")
  openxlsx::writeData(wb, "Module statistics",    as.data.frame(modules$module_stats_df))
  openxlsx::addWorksheet(wb, "Between-module edges")
  openxlsx::writeData(wb, "Between-module edges", as.data.frame(modules$between_module_edges_df))
  openxlsx::addWorksheet(wb, "Module membership")
  openxlsx::writeData(wb, "Module membership",    as.data.frame(modules$membership_df))
  openxlsx::addWorksheet(wb, "All edges")
  openxlsx::writeData(wb, "All edges",            as.data.frame(modules$all_edges_df))
  openxlsx::saveWorkbook(wb, mod_xlsx, overwrite = TRUE)
  message("Saved unweighted modules: ", mod_xlsx)
}

# ---------- Weighted modules ----------
if (run_weighted) {
  message("Running weighted module detection...")
  w_start <- proc.time()
  
  weighted_modules <- detect_PPI_weightedModules(
    gene_list    = dice_genes,
    interactions = phase3_interactions,
    treatment    = opt$treatment,
    control      = opt$control,
    algorithm    = weighted_algo,
    resolution   = weighted_res,
    leiden_itrs  = weighted_leiden_itrs,
    leiden_beta  = weighted_leiden_beta,
    seed         = opt$seed
  )
  
  w_runtime <- proc.time() - w_start
  message("Weighted modules runtime: ", round(w_runtime[["elapsed"]] / 60, 2), " min")
  
  wmod_xlsx <- file.path(opt$outdir, paste0(opt$job_name, "_weighted_modules.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Network modules")
  openxlsx::writeData(wb, "Network modules",       as.data.frame(weighted_modules$network_modules))
  openxlsx::addWorksheet(wb, "Module stats")
  openxlsx::writeData(wb, "Module stats",           as.data.frame(weighted_modules$module_stats))
  openxlsx::addWorksheet(wb, "Jaccard score")
  openxlsx::writeData(wb, "Jaccard score",          as.data.frame(weighted_modules$cmp_jaccard), rowNames = TRUE)
  openxlsx::addWorksheet(wb, "Common genes count")
  openxlsx::writeData(wb, "Common genes count",     as.data.frame(weighted_modules$cmp_count), rowNames = TRUE)
  openxlsx::addWorksheet(wb, "Nodes stats")
  openxlsx::writeData(wb, "Nodes stats",            as.data.frame(weighted_modules$nodes_stats))
  openxlsx::addWorksheet(wb, "Edges")
  openxlsx::writeData(wb, "Edges",                  as.data.frame(weighted_modules$edges_with_modules))
  openxlsx::addWorksheet(wb, "Treatment inter modules")
  openxlsx::writeData(wb, "Treatment inter modules", as.data.frame(weighted_modules$treatment_interMod))
  openxlsx::addWorksheet(wb, "Control inter modules")
  openxlsx::writeData(wb, "Control inter modules",  as.data.frame(weighted_modules$control_interMod))
  openxlsx::saveWorkbook(wb, wmod_xlsx, overwrite = TRUE)
  message("Saved weighted modules: ", wmod_xlsx)
}

# ---------- Total runtime summary ----------
total_runtime <- proc.time() - start_time
total_sec     <- round(total_runtime[["elapsed"]])

message("===================================")
message("Total runtime: ", floor(total_sec / 60), "m ", total_sec %% 60, "s")
message("Finished at: ", Sys.time())
message("===================================")
message("=== Done ===")