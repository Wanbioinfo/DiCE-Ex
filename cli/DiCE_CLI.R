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
  
  library("dplyr")
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
  library("openxlsx")
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
  # NULL -> no rules provided
  if (is.null(x)) return(NULL)
  
  # empty character -> no rules provided
  if (is.character(x) && length(x) == 1) {
    if (is.na(x) || !nzchar(trimws(x))) return(NULL)
    
    x <- tryCatch(
      jsonlite::fromJSON(x, simplifyVector = FALSE),
      error = function(e) die("Failed to parse --dice_rules JSON: ", e$message)
    )
  }
  
  # if config JSON loaded it as data.frame, convert rows to list of rules
  if (is.data.frame(x)) {
    x <- split(x, seq_len(nrow(x)))
    x <- lapply(x, function(row) as.list(row))
  }
  
  # if it is already a list, keep it
  if (!is.list(x) || length(x) == 0) {
    die("dice_rules must be a non-empty list or a valid JSON string.")
  }
  
  x
}

# Map species input to DiCE expected values
normalize_species <- function(x) {
  x <- tolower(trimws(x))
  if (x %in% c("human", "hs", "homo_sapiens", "homo sapiens")) return("human")
  if (x %in% c("mouse", "mm", "mus_musculus", "mus musculus")) return("mouse")
  die("Unsupported species: ", x, " (use 'human' or 'mouse').")
}

# ---------- CLI options ----------
option_list <- list(
  make_option(c("--species"), type="character", default="human", dest="species",
              help="Species: human or mouse [default %default]"),
  
  make_option(c("--dge_file_path"), type="character", default=NULL, dest="dge_file_path",
              help="Path to differential gene expression file"),
  
  make_option(c("--normGeneExp_file_path"), type="character", default=NULL, dest="normGeneExp_file_path",
              help="Path to normalized gene expression matrix"),
  
  make_option(c("--treatment"), type="character", default=NULL, dest="treatment",
              help="Treatment / case group label"),
  
  make_option(c("--control"), type="character", default=NULL, dest="control",
              help="Control / reference group label"),
  
  make_option(c("--loose_criteria"), type="character", default="adj.P.Val", dest="loose_criteria",
              help="Statistical significance metric [default %default]"),
  
  make_option(c("--loose_cutoff"), type="double", default=0.05, dest="loose_cutoff",
              help="Threshold for loose_criteria filtering [default %default]"),
  
  make_option(c("--logFC_cutoff"), type="double", default=0, dest="logFC_cutoff",
              help="Minimum absolute log2 fold-change threshold [default %default]"),
  
  make_option(c("--ig_method"), type="character", default="IG", dest="ig_method",
              help="Information gain method: IG, WIG, none [default %default]"),
  
  make_option(c("--B"), type="integer", default=300, dest="B",
              help="Bootstrap resamples [default %default]"),
  
  make_option(c("--ig_cutoff"), type="character", default="all_mean", dest="ig_cutoff",
              help="IG filtering rule [default %default]"),
  
  make_option(c("--ig_custom_cutoff"), type="double", default=NA, dest="ig_custom_cutoff",
              help="Custom IG threshold when ig_cutoff=custom"),
  
  make_option(c("--corr_method"), type="character", default="pearson", dest="corr_method",
              help="Correlation method [default %default]"),
  
  make_option(c("--centrality_list"), type="character", default="betweenness,eigenvector", dest="centrality_list",
              help="Comma-separated list of centrality metrics"),
  
  make_option(c("--dice_rules"), type = "character",
    default = '[{"type":"centrality","metric":"betweenness","threshold_type":"percent","threshold":25},{"type":"centrality","metric":"eigen vector","threshold_type":"percent","threshold":25}]',
    dest = "dice_rules", help = "JSON string of DiCE selection rules [default %default]"
  ),
  
  make_option(c("--dice_logic"), type="character", default="AND", dest="dice_logic",
              help="How to combine dice_rules: AND or OR [default %default]"),
  
  make_option(c("--run_modules"), action="store_true", default=TRUE, dest="run_modules",
              help="Run PPI module detection after DiCE [default TRUE]"),
  
  make_option(c("--seed"), type="integer", default=123, dest="seed",
              help="Random seed [default %default]"),
  
  make_option(c("--outdir"), type="character", default="dice_cli_output", dest="outdir",
              help="Output directory [default %default]"),
  
  make_option(c("--job_name"), type="character", default="DiCE_CLI_Run", dest="job_name",
              help="Job name prefix for outputs [default %default]"),
  
  make_option(c("--config"), type="character", default=NULL, dest="config",
              help="Optional JSON config file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---------- Load config JSON if provided ----------
if (!is.null(opt$config)) {
  if (!file.exists(opt$config)) die("Config not found: ", opt$config)
  cfg <- jsonlite::fromJSON(opt$config, simplifyVector = FALSE)
  
  message("Config keys: ", paste(names(cfg), collapse = ", "))
  
  for (nm in names(cfg)) opt[[nm]] <- cfg[[nm]]
}

opt$normGeneExp_file_path <- norm_path(opt$normGeneExp_file_path)
opt$dge_file_path         <- norm_path(opt$dge_file_path)

if (!has_file(opt$normGeneExp_file_path)) die("--normGeneExp_file_path missing or file not found.")
if (!has_file(opt$dge_file_path))         die("--dge_file_path missing or file not found.")
if (!has_value(opt$treatment))            die("--treatment is required.")
if (!has_value(opt$control))              die("--control is required.")

species <- normalize_species(opt$species)
set.seed(opt$seed)

ensure_dir(opt$outdir)

# ---------- Load your DiCE code (same as Shiny) ----------
# Adjust this path if your CLI folder is elsewhere.
repo_root <- getwd()
dice_dir  <- file.path(repo_root, "DiCE")

if (!dir.exists(dice_dir)) {
  # fallback: assume you run from repo root
  dice_dir <- "DiCE"
}
if (!dir.exists(dice_dir)) die("Could not find DiCE directory. Expected: ", dice_dir)

# Source all required R files (match your Shiny sources)
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
centrality_list <- trimws(strsplit(opt$centrality_list, ",", fixed = TRUE)[[1]])
centrality_list <- centrality_list[nzchar(centrality_list)]
if (length(centrality_list) == 0) die("No centralities provided in --centrality_list")

# ---------- Run DiCE ----------
message("normGeneExp_file_path: ", opt$normGeneExp_file_path)
message("dge_file_path       : ", opt$dge_file_path)
message("treat/control       : ", opt$treatment, " / ", opt$control)

ig_custom_cutoff <- if (tolower(opt$ig_cutoff) == "custom") opt$ig_custom_cutoff else NULL
if (tolower(opt$ig_cutoff) == "custom" && (is.null(ig_custom_cutoff) || is.na(ig_custom_cutoff))) {
  die("ig_cutoff=custom but --ig_custom_cutoff is missing.")
}

opt$dice_logic <- toupper(trimws(opt$dice_logic))
if (!(opt$dice_logic %in% c("AND", "OR"))) {
  die("--dice_logic must be either 'AND' or 'OR'.")
}

dice_rules <- parse_dice_rules(opt$dice_rules)

# force these (user cannot change)
forced_data_type <- "bulkRNA-seq"
forced_corr_mode <- "directCorr"

dice_res <- perform_DiCE(
  data_type             = "bulkRNA-seq",
  species               = species,
  dge_file_path         = opt$dge_file_path,
  normGeneExp_file_path = opt$normGeneExp_file_path,
  treatment             = opt$treatment,
  control               = opt$control,
  remove_pc_genes       = TRUE,
  loose_criteria        = opt$loose_criteria,
  loose_cutoff          = opt$loose_cutoff,
  logFC_cutoff          = opt$logFC_cutoff,
  ig_method             = opt$ig_method,
  B                     = opt$B,
  ig_cutoff             = opt$ig_cutoff,
  ig_custom_cutoff      = ig_custom_cutoff,
  corr_mode             = "directCorr",
  corr_method           = opt$corr_method,
  centrality_list       = centrality_list,
  dice_rules            = dice_rules,
  dice_logic            = opt$dice_logic
  
)

result_df <- as.data.frame(dice_res$dice_results_df)
phase3_interactions <- as.data.frame(dice_res$interactions_df)

out_xlsx <- file.path(opt$outdir, paste0(opt$job_name, "_dice_results.xlsx"))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "DiCE_results")
openxlsx::writeData(wb, "DiCE_results", result_df)
openxlsx::addWorksheet(wb, "Phase3_interactions")
openxlsx::writeData(wb, "Phase3_interactions", phase3_interactions)
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Saved DiCE results and Phase III interactions: ", out_xlsx)

# ---------- Optional: module detection ----------
if (isTRUE(opt$run_modules)) {
  message("Running unweighted module detection on DiCE genes...")
  
  phase_col <- if ("Phase" %in% names(result_df)) "Phase" else if ("phase" %in% names(result_df)) "phase" else NULL
  if (is.null(phase_col)) die("Could not find Phase column in dice_results_df.")
  
  dice_genes_df <- result_df[result_df[[phase_col]] == "DiCE", , drop = FALSE]
  if (nrow(dice_genes_df) == 0) die("No rows with Phase == 'DiCE' found; cannot run module detection.")
  
  # choose gene symbol column
  gene_col <- c("Gene.Symbol", "Gene.Name", "Gene", "gene")
  gene_col <- gene_col[gene_col %in% names(dice_genes_df)][1]
  if (is.na(gene_col) || is.null(gene_col)) die("Could not find a gene symbol column in results.")
  
  dice_genes <- as.character(dice_genes_df[[gene_col]])
  
  modules <- detect_PPI_unweightedModules(
    gene_list       = dice_genes,
    interactions_df = phase3_interactions,
    seed            = 123
  )
  
  mod_xlsx <- file.path(opt$outdir, paste0(opt$job_name, "_modules.xlsx"))
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "Module summary")
  openxlsx::writeData(wb, "Module summary", as.data.frame(modules$summary_df))
  
  openxlsx::addWorksheet(wb, "Module statistics")
  openxlsx::writeData(wb, "Module statistics", as.data.frame(modules$module_stats_df))
  
  openxlsx::addWorksheet(wb, "Between-module edges")
  openxlsx::writeData(wb, "Between-module edges", as.data.frame(modules$between_module_edges_df))
  
  openxlsx::addWorksheet(wb, "Module membership")
  openxlsx::writeData(wb, "Module membership", as.data.frame(modules$membership_df))
  
  openxlsx::addWorksheet(wb, "All edges")
  openxlsx::writeData(wb, "All edges", as.data.frame(modules$all_edges_df))
  
  openxlsx::saveWorkbook(wb, mod_xlsx, overwrite = TRUE)
  
  message("Saved modules: ", mod_xlsx)
}


message("===================================")
message("Finished at: ", Sys.time())
message("===================================")

message("=== Done ===")