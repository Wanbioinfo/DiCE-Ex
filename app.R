#################################
## app.R – combine all UI tabs ##
#################################

options(shiny.maxRequestSize = 500 * 1024^3)

library(shiny)
library(shinythemes)
library(readr)
library(readxl)
#library(DiCE)
library(future)
library(shinyWidgets)

plan(multisession, workers = 1)
options(future.globals.maxSize = 8 * 1024^3)  # adjust if needed

library(DT)
library(promises)
library(htmltools)
library(visNetwork)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(uuid)
library(jsonlite)
library(zip)
library(vroom)

# Souce DiCE package
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
library('NewWave')
library("praznik")
library("zinbwave")
library("SingleCellExperiment")
library("reticulate")
library("stats")
library("utils")
library("openxlsx")
library("BiocParallel")
library("SummarizedExperiment")
library("parallel")
library("annotate")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("readxl")
library("scuttle")

source("DiCE/annotate_TFs.R")
source("DiCE/calculate_NetCentralities.R")
source("DiCE/calculate_weightedIG.R")
source("DiCE/corr_calculations.R")
source("DiCE/createFinalRanking.R")
source("DiCE/createPPI.R")
source("DiCE/DiCE.R")
source("DiCE/fileReader.R")
source("DiCE/geneExp_denoise_models.R")
source("DiCE/network_analysis.R")
source("DiCE/normalize_df_cols.R")
source("DiCE/Phase_1.R")
source("DiCE/Phase_2.R")
source("DiCE/Phase_3.R")
source("DiCE/Phase_4.R")
source("DiCE/protein_coding_filter.R")
source("DiCE/unweighted_module_analysis.R")
source("DiCE/weighted_module_analysis_helpers.R")
source("DiCE/weighted_module_analysis.R")
source("DiCE/zinbWaVE_denoising.R")
source("DiCE/newWave_denoising.R")

# Source UI pieces
source("ui_home.R")
source("ui_about.R")
source("ui_run.R")
source("ui_results.R")
source("ui_cli.R")
source("ui_team.R")

################################
## Job persistence helpers
################################
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x


jobs_root <- "jobs"
dir.create(jobs_root, showWarnings = FALSE, recursive = TRUE)

job_dir <- function(job_id) file.path(jobs_root, job_id)

write_status <- function(job_id, state, message = "", step = NULL, pct = NULL, job_title = NULL) {
  d <- job_dir(job_id)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  
  jsonlite::write_json(
    list(
      job_id    = job_id,
      job_title = job_title,  
      state     = state,
      message   = message,
      step      = step,
      pct       = pct,
      time      = as.character(Sys.time())
    ),
    path = file.path(d, "status.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

progress_from_dice_msg <- function(msg) {
  if (!is.character(msg) || !nzchar(msg)) return(NULL)
  if (!grepl("Genes in Phase", msg, fixed = TRUE)) return(NULL)
  
  # extract phase number
  phase_num <- suppressWarnings(as.integer(sub(".*Phase\\s*([0-9]+).*", "\\1", msg)))
  if (is.na(phase_num)) return(NULL)
  
  pct_map <- c(`1` = 25, `2` = 50, `3` = 75, `4` = 90)
  pct <- pct_map[as.character(phase_num)]
  if (is.na(pct)) pct <- 60
  
  list(
    pct  = as.numeric(pct),
    step = paste0("Completed DiCE Phase ", phase_num),
    msg  = msg
  )
}

clean_ansi <- function(x) {
  x <- gsub("\033\\[[0-9;]*m", "", x)  # remove ANSI colors like [38;5;232m
  x <- gsub("\\s+", " ", x)
  trimws(x)
}


clean_error_message <- function(x) {
  if (is.null(x)) return(NULL)
  # remove ANSI color codes
  gsub("\\033\\[[0-9;]*m", "", x)
}

friendly_dice_error <- function(err_msg) {
  raw <- clean_ansi(err_msg)
  
  # Default friendly message
  title <- "DiCE could not run with the uploaded files"
  hint  <- "Please check your input files and try again."
  
  # Common patterns (add more as you discover them)
  if (grepl("Gene\\.Name.*not found", raw, ignore.case = TRUE) ||
      grepl("object.*Gene\\.Name.*not found", raw, ignore.case = TRUE)) {
    hint <- paste(
      "Your file does not contain the expected gene identifier column.",
      "Please ensure the DGE file includes a gene column (e.g., Gene.Name / Gene / SYMBOL),",
      "and the expression matrix has genes as columns with sample class label in the last column."
    )
  }
  
  if (grepl("treatment", raw, ignore.case = TRUE) && grepl("control", raw, ignore.case = TRUE)) {
    hint <- paste(
      "The treatment/control labels may not match the class column in your expression matrix.",
      "Please confirm the group labels exactly match the values in the last column."
    )
  }
  
  list(title = title, hint = hint, details = raw)
}

read_status <- function(job_id) {
  f <- file.path(job_dir(job_id), "status.json")
  if (!file.exists(f)) return(NULL)
  jsonlite::read_json(f, simplifyVector = TRUE)
}

job_url <- function(session, job_id) {
  paste0(
    session$clientData$url_protocol, "//",
    session$clientData$url_hostname,
    ifelse(session$clientData$url_port %in% c("", "80", "443"),
           "",
           paste0(":", session$clientData$url_port)),
    session$clientData$url_pathname,
    "?job=", job_id
  )
}

cleanup_jobs <- function(root = "jobs", keep_days = 7) {
  if (!dir.exists(root)) return(invisible(NULL))
  cutoff <- Sys.time() - keep_days * 24 * 3600
  
  job_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  for (d in job_dirs) {
    f <- file.path(d, "status.json")
    if (!file.exists(f)) next
    
    mtime <- file.info(f)$mtime
    if (is.na(mtime) || mtime > cutoff) next
    
    unlink(d, recursive = TRUE, force = TRUE)
  }
}

add_ncbi_gene_links <- function(df, species = "human",
                                gene_col_candidates = c("Gene.Symbol", "Gene.Name", "Gene", "gene")) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  
  # Pick the gene symbol column
  gene_col <- gene_col_candidates[gene_col_candidates %in% names(df)][1]
  if (is.na(gene_col) || is.null(gene_col)) return(df)
  
  symbols <- as.character(df[[gene_col]])
  symbols <- trimws(symbols)
  
  # Choose annotation DB by species
  db <- switch(
    tolower(species),
    "human" = org.Hs.eg.db::org.Hs.eg.db,
    "mouse" = org.Mm.eg.db::org.Mm.eg.db,
    stop("Unsupported species for GeneID mapping: ", species)
  )
  
  # Map SYMBOL -> ENTREZID
  map_df <- AnnotationDbi::select(
    db,
    keys     = unique(symbols),
    keytype  = "SYMBOL",
    columns  = c("ENTREZID")
  )
  
  # If duplicates exist, keep first mapping per SYMBOL
  map_df <- map_df[!duplicated(map_df$SYMBOL), c("SYMBOL", "ENTREZID")]
  id_map <- stats::setNames(map_df$ENTREZID, map_df$SYMBOL)
  
  entrez <- unname(id_map[symbols])
  
  # Build HTML links (fallback to plain symbol if no ID)
  link_prefix <- "https://www.ncbi.nlm.nih.gov/datasets/gene/"
  df[[gene_col]] <- ifelse(
    !is.na(entrez) & nzchar(entrez),
    sprintf('<a href="%s%s/" target="_blank">%s</a>', link_prefix, entrez, symbols),
    symbols
  )
  
  attr(df, "gene_col_linked") <- gene_col
  df
}

################################
## Navbar
################################
make_navbar <- function(selected_tab = "home") {
  navbarPage(
    id    = "main_nav",
    theme = shinytheme("flatly"),
    title = div(
      span("DiCE-Ex", style = "font-weight:600; font-size: 1.3em;"),
      span("v1.2.0",  style = "font-size:0.8em; margin-left:8px; color:#dddddd;")
    ),
    windowTitle = "DiCE-Ex",
    selected    = selected_tab,   # <-- IMPORTANT
    home_tab(),
    about_tab(),
    run_tab(),
    results_tab(),
    cli_tab(), 
    team_tab()
  )
}


################################
## UI
################################

ui <- function(request) {
  
  qs <- parseQueryString(request$QUERY_STRING)
  selected_tab <- if (!is.null(qs$job) && nzchar(qs$job)) "results" else "home"
  
  tagList(
    tags$head(
      tags$style(HTML("
      body { padding-bottom: 50px; }

      .navbar-nav { float: right !important; }
      .navbar-header { float: left !important; }

      .navbar-nav > li > a[data-value='home'] {
        display: none !important;
      }
      
      /* Disable wrapper for downloads */
      .disabled-wrap { pointer-events:none; opacity:0.45; }

      .page-container {
        max-width: 1200px;
        width: 100%;
        margin-left: auto;
        margin-right: auto;
        padding-left: 20px;
        padding-right: 20px;
      }
      
      .home-card-link {
        text-decoration: none !important;
        color: inherit !important;
      }
      
      .feature-card {
        background: #ffffff;
        border-radius: 14px;
        padding: 30px 25px;
        border: 1px solid #e6e6e6;
        text-align: center;         
        transition: all 0.2s ease;
      }
      
      .sample-select-wrap .form-group {
        margin-bottom: 0 !important;
      }
      
      .sample-select-wrap select.form-control {
        height: 51px !important;
        min-height: 51px !important;
        border-radius: 6px !important;
        background-color: #e8f4f1 !important;
        border: 1px solid #b7d7cf !important;
        color: #1f3b4d !important;
        font-size: 15px !important;
        font-weight: 500;
      }
      
      .sample-select-wrap select.form-control:focus {
        border-color: #16a085 !important;
        box-shadow: 0 0 0 0.2rem rgba(22,160,133,0.15) !important;
      }

      .feature-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 6px 18px rgba(0,0,0,0.08);
      }
      
      .feature-icon {
        font-size: 38px;             
        color: #008b7c;
        margin-bottom: 12px;
      }
      
      .feature-card h4 {
        font-weight: 600;
        margin-bottom: 8px;
      }
      
      .feature-card p {
        font-size: 0.95em;
        color: #555;
        margin: 0;
      }



      .navbar .container,
      .navbar .container-fluid {
        max-width: 1200px;
        width: 100%;
        margin-left: auto;
        margin-right: auto;
        padding-left: 20px;
        padding-right: 20px;
      }

      .phase-summary {
        background-color:#f9fafb;
        border-radius:8px;
        border:1px solid #e5e5e5;
        padding:8px 12px;
        margin-top:10px;
        margin-bottom:10px;
        font-size:12px;
        color:#444;
      }
      .phase-summary p {
        margin:0 0 2px 0;
        font-family:monospace;
      }

      /* --- PPI subnetwork panel --- */
      .module-network-panel {
        background-color: #f8fafc;
        border-radius: 10px;
        border: 1px solid #d7d8db;
        padding: 10px;
        position: relative;
        min-height: 300px;
      }

      #ppi_module_network {
        width: 100% !important;
        position: relative;
      }
      #ppi_module_network .vis-network {
        width: 100% !important;
      }

      #ppi_module_network .vis-network-export {
        position: absolute !important;
        top: 10px;
        right: 10px;
        z-index: 20;
      }

      .modal-content { border-radius: 12px; }
      .modal-body { font-size: 16px; text-align: center; }

      @media (max-width: 767px) {
        .navbar-header { float: none !important; }
        .navbar-nav    { float: none !important; }
      }
    ")),
      tags$script(HTML("
      // click on brand -> go home
      $(document).on('click', '.navbar-brand', function(e){
        e.preventDefault();
        Shiny.setInputValue('brand_click', new Date().getTime());
      });

      // smooth scroll to docs section in About tab
      Shiny.addCustomMessageHandler('scroll_to_docs', function(message) {
        setTimeout(function() {
          var el = document.getElementById('docs_section');
          if (el) el.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }, 300);
      });

      Shiny.addCustomMessageHandler('scroll_to_updates', function(message) {
        setTimeout(function() {
          var el = document.getElementById('updates_section');
          if (el) el.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }, 300);
      });

      // visually show sample filenames next to fileInputs
      Shiny.addCustomMessageHandler('loadFiles', function(msg) {
        if (msg.expr) $('#expr_file').parent().find('.form-control').val(msg.expr);
        if (msg.dge)  $('#dge_file').parent().find('.form-control').val(msg.dge);
      });

      // external trigger to call visNetwork exportPNG()
      Shiny.addCustomMessageHandler('trigger_vis_export', function(message) {
        if (window.network) window.network.exportPNG();
      });
  
    "))
    ),
    
    make_navbar(selected_tab),   
    
    tags$footer(
      class = "app-footer",
      style = "
      position: fixed;
      left: 0; right: 0; bottom: 0;
      padding: 6px 20px;
      background-color: #003a5d;
      color: #ffffff;
      font-size: 0.9em;
    ",
      div(
        style = "text-align: center; line-height: 1.4;",
        span("@ Wan Lab, Indiana University School of Medicine · 2025"),
        tags$br(),
        span(
          "This website is free and open to all users and there is no login requirement.",
          style = "font-size:1.1em;"
        ),
        span("  ·  "),
        a(
          "Lab Website",
          href   = "https://wanbioinfo.github.io/Lab/",
          target = "_blank",
          style  = "font-size:0.9em; color:#ffffff; text-decoration:underline;"
        ),
        span("  ·  "),
        a(
          "License",
          href   = "https://opensource.org/licenses/MIT",
          target = "_blank",
          style  = "font-size:0.9em; color:#ffffff; text-decoration:underline;"
        )
      )
    )
  )
}



################################
## Server
################################
server <- function(input, output, session) {
  
  # Clean up old jobs once per session start
  cleanup_jobs(keep_days = 1)
  
  is_phase_message <- function(x) {
    is.character(x) && length(x) == 1 && grepl("^#Genes in Phase", x)
  }
  
  job_id_qs <- reactive({
    qs <- parseQueryString(session$clientData$url_search)
    if (is.null(qs$job) || !nzchar(qs$job)) return(NULL)
    qs$job
  })
  
  observeEvent(job_id_qs(), {
    jid <- job_id_qs()
    if (is.null(jid)) return()
    updateNavbarPage(session, "main_nav", "results")
  }, ignoreInit = FALSE)
  
  job_status <- reactivePoll(
    intervalMillis = 2000,
    session = session,
    checkFunc = function() {
      jid <- job_id_qs()
      if (is.null(jid)) return(NA_real_)
      f <- file.path(job_dir(jid), "status.json")
      if (!file.exists(f)) return(NA_real_)
      as.numeric(file.info(f)$mtime)
    },
    valueFunc = function() {
      jid <- job_id_qs()
      if (is.null(jid)) return(NULL)
      read_status(jid)
    }
  )

  
  ################################
  ## 1. Navigation
  ################################
  observeEvent(input$brand_click, {
    updateNavbarPage(session, "main_nav", "home")
  })
  observeEvent(input$go_run, {
    updateNavbarPage(session, "main_nav", "run")
  })
  observeEvent(input$home_docs, {
    updateNavbarPage(session, "main_nav", "about")
    session$sendCustomMessage("scroll_to_docs", list())
  })
  observeEvent(input$home_updates, {
    updateNavbarPage(session, "main_nav", "about")
    session$sendCustomMessage("scroll_to_updates", list())
  })
  observeEvent(input$go_about_tab, {
    updateNavbarPage(session, "main_nav", selected = "about")
  })
  
  
  output$expr_loaded_badge <- renderUI({
    req(expr_df())
    tags$div(
      style="margin-top:6px; color:#198754; font-size:13px;",
      icon("check-circle"), " Sample expression data loaded"
    )
  })
  
  output$dge_loaded_badge <- renderUI({
    req(dge_df())
    tags$div(
      style="margin-top:6px; color:#198754; font-size:13px;",
      icon("check-circle"), " Sample DGE data loaded"
    )
  })
  
  
  ################################
  ## 2. Log system
  ################################
  dice_log <- reactiveVal("")
  add_log <- function(msg) {
    old <- isolate(dice_log())
    msg <- paste(msg, collapse = "\n")
    new <- paste(old, msg, sep = ifelse(old == "", "", "\n"))
    dice_log(new)
  }
  output$dice_log <- renderText(dice_log())
  current_status <- reactiveVal("Waiting…")
  output$dice_modal_status <- renderText(current_status())

  ################################
  ## 3. Results + state
  ################################
  dice_result   <- reactiveVal(NULL)
  phase_lines   <- reactiveVal(character(0))
  job_title_val <- reactiveVal(NULL)
  
  modules_summary    <- reactiveVal(NULL)
  modules_membership <- reactiveVal(NULL)
  modules_edges      <- reactiveVal(NULL)
  
  ################################
  ## 3b. Download readiness flag  <-- PUT IT HERE
  ################################
  downloads_ready <- reactive({
    !is.null(dice_result()) &&
      is.data.frame(dice_result()) &&
      nrow(dice_result()) > 0
  })
  
  modules_ready <- reactive({
    !is.null(modules_summary()) &&
      is.data.frame(modules_summary()) &&
      nrow(modules_summary()) > 0 &&
      !is.null(modules_membership()) &&
      is.data.frame(modules_membership()) &&
      nrow(modules_membership()) > 0
  })
  
  output$download_buttons_ui <- renderUI({
    if (isTRUE(downloads_ready())) {
      downloadButton("download_dice_results", "Download DiCE results", class = "btn btn-success")
    } else {
      div(
        class = "disabled-wrap",
        downloadButton("download_dice_results", "Download DiCE results", class = "btn btn-success"),
        tags$div(style="margin-top:6px; font-size:12px; color:#666;",
                 "Downloads will be enabled after results finish loading.")
      )
    }
  })
  
  output$download_modules_ui <- renderUI({
    if (isTRUE(modules_ready())) {
      downloadButton("download_modules_both", "Download modules", class = "btn btn-success")
    } else {
      div(
        class = "disabled-wrap",
        downloadButton("download_modules_both", "Download modules", class = "btn btn-success"),
        tags$div(style="margin-top:6px; font-size:12px; color:#666;",
                 "Module downloads will be enabled after module results finish loading.")
      )
    }
  })
  
  
  expr_df <- reactiveVal(NULL)
  dge_df  <- reactiveVal(NULL)
  

  
  ################################
  ## 4. Job link + status polling (NO reload)
  ################################
  
  output$job_status_ui <- renderUI({
    st <- job_status()
    if (is.null(st)) return(NULL)
    
    msg <- st$message
    show_msg <- is.character(msg) &&
      length(msg) == 1 &&
      !is.na(msg) &&
      nzchar(msg) &&
      !is_phase_message(msg)   # 👈 FILTER HERE
    
    div(
      class = "phase-summary",
      tags$p(paste("Job:", st$job_id %||% "")),
      tags$p(paste("Status:", st$state %||% "unknown")),
      if (show_msg) tags$p(msg)
    )
  })
  
  output$error_help_ui <- renderUI({
  st <- job_status()
  if (is.null(st) || st$state != "error") return(NULL)

  div(
    class = "phase-summary",
    style = "border-left:5px solid #dc3545; background:#fff5f5;",

    tags$h4(
      icon("exclamation-triangle"),
      " DiCE failed — please check your input data",
      style = "color:#b02a37;"
    ),

    if (!is.null(st$message) && nzchar(st$message)) {
      tags$p(
        tags$b("Error message:"),
        tags$br(),
        tags$code(st$message)
      )
    },

    tags$hr(),

    tags$p(tags$b("Common input issues to check:")),

    tags$ul(
      tags$li(
        tags$b("Expression matrix: "),
        "Genes must be columns, samples must be rows, and the class/phenotype label must be in the last column.",
        "Do not place sample names as the first column."
      ),
      tags$li(
        tags$b("DGE file: "),
        "Must contain a valid gene identifier column (e.g., Gene or Gene.Name), ",
        "logFC, p-value, and adjusted p-value."
      ),
      tags$li(
        tags$b("Gene identifiers: "),
        "Must match between expression and DGE files."
      ),
      tags$li(
        tags$b("Group labels: "),
        "Must exactly match the values in the class column (case-sensitive)."
      )
    )
  )
})

  
  output$job_progress_ui <- renderUI({
    st <- job_status()
    if (is.null(st)) return(NULL)
    
    # ---- robust coercion (safe for NULL / NA / length-0) ----
    state <- st$state
    if (is.null(state) || length(state) == 0) state <- "unknown"
    state <- as.character(state)[1]
    if (is.na(state) || !nzchar(state)) state <- "unknown"
    
    # hide progress box once job is finished or errored
    if (identical(state, "finished") || identical(state, "error")) {
      return(NULL)
    }
    
    step <- st$step
    if (is.null(step) || length(step) == 0) step <- "Waiting…"
    step <- as.character(step)[1]
    if (is.na(step) || !nzchar(step)) step <- "Waiting…"
    
    msg <- st$message
    if (is.null(msg) || length(msg) == 0) msg <- ""
    msg <- as.character(msg)[1]
    
    show_msg <- !is.na(msg) &&
      nzchar(msg) &&
      !is_phase_message(msg)  
    
    
    pct <- suppressWarnings(as.numeric(st$pct))
    if (is.null(pct) || length(pct) == 0 || is.na(pct)) pct <- 0
    pct <- max(0, min(100, pct))
    # ---------------------------------------------------------
    
    div(
      class = "phase-summary",
      tags$p(tags$b("Progress")),
      tags$p(paste0("Status: ", state)),
      tags$p(paste0("Step: ", step)),
      tags$div(
        style = "background:#e9ecef; border-radius:10px; overflow:hidden; height:14px; margin-top:6px;",
        tags$div(
          style = paste0("width:", pct, "%; height:14px; background:#0d6efd;")
        )
      ),
      tags$p(style = "margin-top:6px;", paste0("Completed: ", pct, "%")),
      if (show_msg) tags$p(msg)
    )
  })
  
  
  loaded_once <- reactiveVal(FALSE)
  
  observeEvent(job_status(), {
    st  <- job_status()
    jid <- job_id_qs()
    if (is.null(st) || is.null(jid)) return()
    
    if (!is.null(st$job_title) &&
        is.character(st$job_title) &&
        length(st$job_title) == 1 &&
        !is.na(st$job_title) &&
        nzchar(st$job_title)) {
      job_title_val(st$job_title)
    }
    
    updateNavbarPage(session, "main_nav", "results")
  
    # Load outputs once when finished
    if (identical(st$state, "finished") && !isTRUE(loaded_once())) {
      jd <- job_dir(jid)
      
      if (file.exists(file.path(jd, "dice_result.rds")))
        dice_result(readRDS(file.path(jd, "dice_result.rds")))
      if (file.exists(file.path(jd, "modules_summary.rds")))
        modules_summary(readRDS(file.path(jd, "modules_summary.rds")))
      if (file.exists(file.path(jd, "modules_membership.rds")))
        modules_membership(readRDS(file.path(jd, "modules_membership.rds")))
      if (file.exists(file.path(jd, "modules_edges.rds")))
        modules_edges(readRDS(file.path(jd, "modules_edges.rds")))
      if (file.exists(file.path(jd, "expr_input.rds")))
        expr_df(readRDS(file.path(jd, "expr_input.rds")))
      if (file.exists(file.path(jd, "dge_input.rds")))
        dge_df(readRDS(file.path(jd, "dge_input.rds")))
      
      
      log_file <- file.path(jd, "log.txt")
      if (file.exists(log_file)) {
        log_lines <- readLines(log_file, warn = FALSE)
        phase_lines(grep("Genes in Phase", log_lines, value = TRUE))
        dice_log("")
        add_log("----- Loaded saved job log -----")
        add_log(log_lines)
      }
      
      loaded_once(TRUE)
      add_log("Loaded finished results from bookmarked job link.")
    }
  }, ignoreInit = FALSE)
  
  ################################
  ## 5. Upload paths (uploads OR sample files)
  ################################
  expr_path <- reactiveVal(NULL)
  dge_path  <- reactiveVal(NULL)
  
  observeEvent(input$expr_file, {
    req(input$expr_file)
    path <- input$expr_file$datapath
    expr_path(path)
    
    ext <- tolower(tools::file_ext(input$expr_file$name))
    expr_tbl <- switch(
      ext,
      "csv"  = readr::read_csv(path, show_col_types = FALSE),
      "tsv"  = readr::read_tsv(path,  show_col_types = FALSE),
      "txt"  = readr::read_tsv(path,  show_col_types = FALSE),
      "xls"  = readxl::read_excel(path),
      "xlsx" = readxl::read_excel(path),
      "rds"  = {
        obj <- readRDS(path)
        if (is.matrix(obj)) obj <- as.data.frame(obj)
        if (!is.data.frame(obj)) stop("Expression RDS must contain a data frame or matrix.")
        obj
      },
      readr::read_delim(path, delim = NULL, show_col_types = FALSE)
    )
    expr_df(as.data.frame(expr_tbl))
  })
  
  observeEvent(input$dge_file, {
    req(input$dge_file)
    path <- input$dge_file$datapath
    dge_path(path)
    
    ext <- tolower(tools::file_ext(input$dge_file$name))
    dge_tbl <- switch(
      ext,
      "csv"  = readr::read_csv(path, show_col_types = FALSE),
      "tsv"  = readr::read_tsv(path,  show_col_types = FALSE),
      "txt"  = readr::read_tsv(path,  show_col_types = FALSE),
      "xls"  = readxl::read_excel(path),
      "xlsx" = readxl::read_excel(path),
      "rds"  = {
        obj <- readRDS(path)
        if (is.matrix(obj)) obj <- as.data.frame(obj)
        if (!is.data.frame(obj)) stop("DGE RDS must contain a data frame or matrix.")
        obj
      },
      readr::read_delim(path, delim = NULL, show_col_types = FALSE)
    )
    dge_df(as.data.frame(dge_tbl))
  })
  
  ################################
  ## 6. Update phase selector
  ################################
  observeEvent(dice_result(), {
    df <- dice_result()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    
    phase_col <- NULL
    if ("Phase" %in% names(df)) phase_col <- "Phase"
    if ("phase" %in% names(df)) phase_col <- "phase"
    
    if (!is.null(phase_col)) {
      phases <- sort(unique(as.character(df[[phase_col]])))
      updateSelectInput(
        session,
        "dice_phase_filter",
        choices  = c("All phases" = "all", stats::setNames(phases, phases)),
        selected = "all"
      )
    }
  })
  
  ################################
  ## 7. Results UI outputs
  ################################
  output$dice_log_ui <- renderText({
    log_txt <- dice_log()
    req(log_txt)
    log_txt
  })
  
  output$dice_job_title <- renderUI({
    jt <- job_title_val()
    if (is.null(jt) || jt == "") return(NULL)
    tags$h4(paste("Job:", jt), style = "margin-top:5px; margin-bottom:5px; color:#333;")
  })
  
  output$phase_summary <- renderUI({
    pl <- phase_lines()
    if (!length(pl)) return(NULL)
    tags$div(class = "phase-summary", lapply(pl, function(x) tags$p(x)))
  })
  
  ################################
  ## 8. DiCE results table
  ################################
  dice_table_data <- reactive({
    df <- dice_result()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    
    phase_col <- NULL
    if ("Phase" %in% names(df)) phase_col <- "Phase"
    if ("phase" %in% names(df)) phase_col <- "phase"
    
    if (!is.null(input$dice_phase_filter) &&
        input$dice_phase_filter != "all" &&
        !is.null(phase_col)) {
      df <- df[df[[phase_col]] == input$dice_phase_filter, , drop = FALSE]
    }
    
    q <- input$dice_search
    q <- if (is.null(q) || length(q) == 0 || is.na(q)) "" else trimws(as.character(q))
    
    if (nzchar(q)) {
      gene_col <- dplyr::case_when(
        "Gene.Name" %in% names(df) ~ "Gene.Name",
        "Gene"      %in% names(df) ~ "Gene",
        "gene"      %in% names(df) ~ "gene",
        TRUE ~ NA_character_
      )
      if (!is.na(gene_col)) {
        keep <- grepl(q, as.character(df[[gene_col]]), ignore.case = TRUE)
        df   <- df[keep, , drop = FALSE]
      }
    }
    # Now add hyperlinks
    df <- add_ncbi_gene_links(df, species = tolower(input$species))
    
    df
  })
  
  output$dice_table <- DT::renderDT({
    df <- dice_table_data()
    req(df)
    
    # Identify which column contains HTML links
    gene_col <- attr(df, "gene_col_linked")
    gene_idx <- if (!is.null(gene_col) && gene_col %in% names(df)) which(names(df) == gene_col) else NULL
    
    # Escape everything except the gene column
    escape_cols <- if (!is.null(gene_idx)) setdiff(seq_along(df), gene_idx) else TRUE
    
    DT::datatable(
      df,
      filter   = "none",
      rownames = FALSE,
      escape   = escape_cols,
      options  = list(dom = "lrtip", pageLength = 20, scrollX = TRUE),
      selection = "none"
    )
  })
  
  ################################
  ## 9. Module tables
  ################################
  output$modules_summary_table <- DT::renderDT({
    df <- modules_summary()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    DT::datatable(df, filter="none", rownames=FALSE,
                  options=list(scrollX=TRUE, pageLength=10, dom="lrtip"))
  })
  
  observeEvent(modules_membership(), {
    df <- modules_membership()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    mods <- sort(unique(df$Module))
    updateSelectInput(session, "mm_module_filter",
                      choices=c("All modules"="all", setNames(as.character(mods), mods)),
                      selected="all")
  })
  
  membership_filtered <- reactive({
    df <- modules_membership()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    
    if (!is.null(input$mm_module_filter) && input$mm_module_filter != "all") {
      df <- df[df$Module == as.numeric(input$mm_module_filter), , drop = FALSE]
    }
    
    q <- input$mm_gene_search
    q <- if (is.null(q) || length(q) == 0 || is.na(q)) "" else trimws(as.character(q))
    
    if (nzchar(q)) {
      df <- df[grepl(q, df$Gene, ignore.case = TRUE), , drop = FALSE]
    }
    df
  })
  
  
  output$modules_membership_table <- DT::renderDT({
    df <- membership_filtered()
    req(df)
    DT::datatable(
      df,
      filter="none",
      rownames=FALSE,
      selection=list(mode="single", target="row"),
      options=list(scrollX=TRUE, pageLength=20, dom="lrtip")
    )
  })
  
  mem_proxy <- DT::dataTableProxy("modules_membership_table")
  
  observeEvent(input$mm_gene_search, {
    q <- input$mm_gene_search
    q <- if (is.null(q) || length(q) == 0 || is.na(q)) "" else trimws(as.character(q))
    
    if (!nzchar(q)) {
      DT::selectRows(mem_proxy, NULL)
    }
  })
  
  
  ################################
  ## 10. Downloads
  ################################
  output$download_dice_results <- downloadHandler(
    filename = function() {
      job <- gsub("[^A-Za-z0-9_-]", "_", input$job_title)
      paste0("DiCE_results_", job, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".xlsx")
    },
    content = function(file) {
      
      df <- dice_result()   # ORIGINAL clean dataframe
      req(df)
      
      if (!is.data.frame(df)) df <- as.data.frame(df)
      
      openxlsx::write.xlsx(df, file, rowNames = FALSE)
    }
  )
  
  output$download_modules_both <- downloadHandler(
    filename = function() {
      job <- gsub("[^A-Za-z0-9_-]", "_", input$job_title)
      paste0("DiCE_modules_", job, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".xlsx")
    },
    content = function(file) {
      s <- modules_summary()
      m <- modules_membership()
      req(s, m)
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Module summary")
      openxlsx::writeData(wb, "Module summary", as.data.frame(s))
      openxlsx::addWorksheet(wb, "Module membership")
      openxlsx::writeData(wb, "Module membership", as.data.frame(m))
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$download_gene_plot <- downloadHandler(
    filename = function() paste0(input$ppi_clicked_gene, "_expression.png"),
    content = function(file) {
      req(input$ppi_clicked_gene)
      g   <- input$ppi_clicked_gene
      dat <- expr_df()
      cond_col <- names(dat)[ncol(dat)]
      
      df_long <- data.frame(
        Condition  = dat[[cond_col]],
        Expression = dat[[g]]
      )
      
      p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Condition, y = Expression, fill = Condition)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
        ggplot2::labs(title = paste("Expression of", g), x = "", y = "Normalized expression") +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::theme(legend.position = "none")
      
      ggsave(file, plot = p, width = 7, height = 5, dpi = 300)
    }
  )
  
  ###########################################
  ## 11. Network
  ###########################################
  selected_module_id <- reactive({
    if (is.null(input$mm_module_filter) || input$mm_module_filter == "all") return("all")
    as.numeric(input$mm_module_filter)
  })
  
  selected_gene_id <- reactive({
    mem <- membership_filtered()
    sel <- input$modules_membership_table_rows_selected
    if (length(sel) == 1) mem$Gene[sel] else NULL
  })
  
  module_nodes <- reactive({
    mem <- modules_membership()
    req(mem)
    if (!is.data.frame(mem)) mem <- as.data.frame(mem)
    mid <- selected_module_id()
    if (identical(mid, "all")) return(mem)
    mem[mem$Module == mid, , drop = FALSE]
  })
  
  module_edges <- reactive({
    edges_list <- modules_edges()
    req(edges_list)
    mid <- selected_module_id()
    
    if (identical(mid, "all")) {
      all_edges <- bind_rows(lapply(edges_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
      if (nrow(all_edges) == 0) return(data.frame(Gene1=character(), Gene2=character(), stringsAsFactors=FALSE))
      return(all_edges)
    }
    
    ed <- edges_list[[as.character(mid)]]
    if (is.null(ed)) return(data.frame(Gene1=character(), Gene2=character(), stringsAsFactors=FALSE))
    as.data.frame(ed, stringsAsFactors = FALSE)
  })
  
  output$ppi_module_network <- renderVisNetwork({
    nodes <- module_nodes()
    edges <- module_edges()
    req(nodes, edges)
    
    nodes_df <- data.frame(
      id    = nodes$Gene,
      label = nodes$Gene,
      title = paste0(nodes$Gene, "<br>Module: ", nodes$Module, "<br>Degree: ", nodes$Degree_inModule),
      value = nodes$Degree_inModule,
      stringsAsFactors = FALSE
    )
    edges_df <- data.frame(from = edges$Gene1, to = edges$Gene2, stringsAsFactors = FALSE)
    
    if (nrow(nodes_df) == 0 || nrow(edges_df) == 0) {
      return(visNetwork(data.frame(), data.frame()) %>% visOptions(nodesIdSelection = FALSE))
    }
    
    highlight_gene <- selected_gene_id()
    
    if (!is.null(highlight_gene) && nzchar(highlight_gene)) {
      keep_edges <- edges_df$from == highlight_gene | edges_df$to == highlight_gene
      edges_sub  <- edges_df[keep_edges, , drop = FALSE]
      if (nrow(edges_sub) > 0) {
        keep_nodes <- unique(c(edges_sub$from, edges_sub$to))
        nodes_df <- nodes_df[nodes_df$id %in% keep_nodes, , drop = FALSE]
        edges_df <- edges_sub
      }
      
      nodes_df$color.background <- ifelse(nodes_df$id == highlight_gene, "#e74c3c", "#95a5a6")
      nodes_df$color.border     <- ifelse(nodes_df$id == highlight_gene, "#c0392b", "#7f8c8d")
    }
    
    max_nodes <- 200
    if (is.null(highlight_gene) && nrow(nodes_df) > max_nodes) {
      top_ids  <- nodes_df$id[order(-nodes_df$value)][1:max_nodes]
      nodes_df <- nodes_df[nodes_df$id %in% top_ids, , drop = FALSE]
      edges_df <- edges_df[edges_df$from %in% nodes_df$id & edges_df$to %in% nodes_df$id, , drop = FALSE]
    }
    
    file_name <- if (!is.null(highlight_gene) && nzchar(highlight_gene)) {
      paste0(highlight_gene, "_module_network")
    } else {
      paste0("Module_", selected_module_id(), "_network")
    }
    
    visNetwork(nodes_df, edges_df) %>%
      visNodes(shape = "dot", size = 32, font = list(size = 28)) %>%
      visOptions(highlightNearest = FALSE, nodesIdSelection = FALSE) %>%
      visIgraphLayout(layout = "layout_with_fr") %>%
      visPhysics(enabled = FALSE) %>%
      visEvents(
        selectNode = "function(nodes) {
          if (nodes.nodes.length > 0) {
            Shiny.setInputValue('ppi_clicked_gene', nodes.nodes[0], {priority: 'event'});
          }
        }",
        afterDrawing = "function() { window.network = this.network; }"
      ) %>%
      visExport(type = "png", name = file_name)
  })
  
  observeEvent(input$export_network_png, {
    session$sendCustomMessage("trigger_vis_export", list())
  })
  
  ################################
  ## 12. Download sample data (.zip)
  ################################
  output$download_sample_zip <- downloadHandler(
    filename = function() "sample_data.zip",
    content = function(file) {
      zip::zip(
        zipfile = file,
        files = c("sample_data/sample_Human_data_geneExp.csv",
                  "sample_data/sample_Human_data_DGE.csv",
                  "sample_data/sample_Mouse_data_geneExp.csv",
                  "sample_data/sample_Mouse_data_DGE.csv")
      )
    }
  )
  
  ################################
  ## 13. Load sample data
  ################################
  ################################
  ## 13. Load sample data
  ################################
  observeEvent(input$sample_dataset, {
    
    req(input$sample_dataset)
    if (input$sample_dataset == "") return()
    
    if (input$sample_dataset == "human") {
      
      updateTextInput(session, "job_title",     value = "Human_DiCE_test")
      updateTextInput(session, "group_treat",   value = "Tumor")
      updateTextInput(session, "group_control", value = "Normal")
      updateRadioButtons(session, "species", selected = "Human")
      
      sample_expr <- "sample_data/sample_Human_data_geneExp.csv"
      sample_dge  <- "sample_data/sample_Human_data_DGE.csv"
      
      expr_path(sample_expr)
      dge_path(sample_dge)
      
      expr_df(vroom::vroom(sample_expr, show_col_types = FALSE))
      dge_df(vroom::vroom(sample_dge, show_col_types = FALSE))
      
      session$sendCustomMessage(
        "loadFiles",
        list(
          expr = "sample_Human_data_geneExp.csv",
          dge  = "sample_Human_data_DGE.csv"
        )
      )
      
      showModal(modalDialog(
        title = "Human sample data loaded",
        HTML(paste0(
          "Human sample expression and DGE files have been loaded successfully.<br><br>",
          "<b>Job title:</b> Human_DiCE_test<br>",
          "<b>Treatment:</b> Tumor<br>",
          "<b>Control:</b> Normal<br>",
          "<b>Species:</b> Human"
        )),
        easyClose = TRUE,
        footer = modalButton("Dismiss")
      ))
      
    } else if (input$sample_dataset == "mouse") {
      
      updateTextInput(session, "job_title",     value = "Mouse_DiCE_test")
      updateTextInput(session, "group_treat",   value = "POR")
      updateTextInput(session, "group_control", value = "WT")
      updateRadioButtons(session, "species", selected = "Mouse")
      
      sample_expr <- "sample_data/sample_Mouse_data_geneExp.csv"
      sample_dge  <- "sample_data/sample_Mouse_data_DGE.csv"
      
      expr_path(sample_expr)
      dge_path(sample_dge)
      
      expr_df(vroom::vroom(sample_expr, show_col_types = FALSE))
      dge_df(vroom::vroom(sample_dge, show_col_types = FALSE))
      
      session$sendCustomMessage(
        "loadFiles",
        list(
          expr = "sample_Mouse_data_geneExp.csv",
          dge  = "sample_Mouse_data_DGE.csv"
        )
      )
      
      showModal(modalDialog(
        title = "Mouse sample data loaded",
        HTML(paste0(
          "Mouse sample expression and DGE files have been loaded successfully.<br><br>",
          "<b>Job title:</b> Mouse_DiCE_test<br>",
          "<b>Treatment:</b> POR<br>",
          "<b>Control:</b> WT<br>",
          "<b>Species:</b> Mouse"
        )),
        easyClose = TRUE,
        footer = modalButton("Dismiss")
      ))
    }
  })
  
  ################################
  ## 14. Run DiCE (creates bookmark link)
  ################################
  observeEvent(input$run_dice, {
    
    missing <- character()
    
    if (is.null(input$job_title) || trimws(input$job_title) == "")
      missing <- c(missing, "Job title")
    if (is.null(expr_path()))
      missing <- c(missing, "Gene expression data file")
    if (is.null(dge_path()))
      missing <- c(missing, "Differential gene expression results file")
    if (is.null(input$group_treat) || trimws(input$group_treat) == "")
      missing <- c(missing, "Treatment / case group label")
    if (is.null(input$group_control) || trimws(input$group_control) == "")
      missing <- c(missing, "Control / reference group label")
    
    if (length(missing) > 0) {
      showModal(modalDialog(
        title = "Missing required input",
        easyClose = TRUE,
        footer = modalButton("Close"),
        HTML(paste0(
          "Please fill in the following before running DiCE:<br><br>",
          paste("&bull; ", missing, collapse = "<br>")
        ))
      ))
      return()
    }
    
    # Create job id + persist inputs
    job_id <- uuid::UUIDgenerate()
    jd <- job_dir(job_id)
    dir.create(jd, showWarnings = FALSE, recursive = TRUE)
    
    saveRDS(expr_df(), file.path(jd, "expr_input.rds"))
    saveRDS(dge_df(),  file.path(jd, "dge_input.rds"))
    
    write_status(job_id, "queued", "Job submitted.", job_title = input$job_title)
    
    link <- job_url(session, job_id)
    
    # Reset UI state
    dice_log("")
    phase_lines(character(0))
    job_title_val(input$job_title)
    add_log(paste0("Starting DiCE run. Job link: ", link))
    
    params <- list(
      job_title    = input$job_title,
      species      = tolower(input$species),
      treat        = input$group_treat,
      control      = input$group_control,
      sig_metric   = input$significant_metric,
      sig_thresh   = input$significant_thresh,
      logfc_thresh = input$phase1_logfc_thresh,
      ig_method    = input$phase2_ig_method,
      B            = input$phase2_B,
      ig_cutoff    = input$phase2_ig_cutoff,
      custom_ig    = input$phase2_custom_ig,
      corr_type    = tolower(input$phase3_corr),
      centralities = input$phase4_cents,
      dice_cutoff  = if (input$phase4_dice_cutoff == "topK")
        paste0("top", input$phase4_topK) else input$phase4_dice_cutoff
    )
    
    dge_file  <- dge_path()
    expr_file <- expr_path()
    
    current_status("Initializing DiCE run…")
    
    # One modal only (includes link)
    showModal(modalDialog(
      easyClose = FALSE,
      footer    = NULL,
      size      = "m",
      div(
        style = "text-align:center;",
        icon("spinner", class = "fa-spin fa-3x", style = "margin-bottom:10px;"),
        h4("Running DiCE… Please wait."),
        HTML(paste0(
          "<div style='margin-top:12px; text-align:left;'>",
          "<b>Bookmarkable results link:</b><br>",
          "<a href='", link, "' target='_blank'>", link, "</a><br>",
          "<span style='font-size:12px; color:#555;'>This page will auto-update.</span>",
          "</div>"
        ))
      )
    ))
    
    future::future({
      
      write_status(job_id, "running",
                   message="DiCE is running...",
                   step="Running DiCE (Phases I–IV)…",
                   pct=10,
                   job_title=params$job_title)
      
      jd <- job_dir(job_id)
      
      log_vec   <- character()
      result_df <- NULL
      modules   <- NULL
      dice_res <- NULL
      phase3_interactions <- NULL
      
      log_vec <- capture.output(
        withCallingHandlers(
          {
            dice_res <- perform_DiCE(
              data_type             = "bulkRNA-seq",
              species               = params$species,
              dge_file_path         = dge_file,
              normGeneExp_file_path = expr_file,
              treatment             = params$treat,
              control               = params$control,
              loose_criteria        = params$sig_metric,
              loose_cutoff          = params$sig_thresh,
              logFC_cutoff          = params$logfc_thresh,
              ig_method             = params$ig_method,
              B                     = params$B,
              ig_cutoff             = params$ig_cutoff,
              ig_custom_cutoff      = params$custom_ig,
              corr_mode             = "directCorr",
              corr_method           = params$corr_type,
              centrality_list       = params$centralities,
              min_passCount         = length(params$centralities),
              cutoff                = params$dice_cutoff
            )
            
            
            write_status(job_id, "running",
                         message="DiCE finished, running module detection…",
                         step="Detecting PPI modules…",
                         pct=95,
                         job_title=params$job_title)
            
            result_df <- as.data.frame(dice_res$dice_results_df)
            
            phase3_interactions <- dice_res$interactions_df
            
            phase_col <- if ("Phase" %in% names(result_df)) "Phase" else if ("phase" %in% names(result_df)) "phase" else NULL
            if (is.null(phase_col)) stop("Could not find Phase column in DiCE results.")
            
            dice_genes_df <- result_df[result_df[[phase_col]] == "DiCE", , drop = FALSE]
            if (nrow(dice_genes_df) == 0) stop("No rows with Phase == 'DiCE' found; cannot run module detection.")
            
            print(class(dice_genes_df))
            print(colnames(dice_genes_df))
            dice_genes <- dice_genes_df$Gene.Symbol
            
            modules <- detect_PPI_unweightedModules(
              gene_list     = dice_genes,
              species       = params$species,
              seed          = 123
            )
            write_status(job_id, "running",
                         message="Saving outputs…",
                         step="Saving outputs…",
                         pct=98,
                         job_title=params$job_title)
            
            write_status(job_id, "finished",
                         message="Done.",
                         step="Finished.",
                         pct=100,
                         job_title=params$job_title)
            
            
          },

          message = function(m) {
            msg <- conditionMessage(m)
            cat(msg, "\n")  # keep your log
            
            # LIVE progress update when DiCE prints "Genes in Phase ..."
            p <- progress_from_dice_msg(msg)
            if (!is.null(p)) {
              write_status(
                job_id    = job_id,
                state     = "running",
                message   = p$msg,
                step      = p$step,
                pct       = p$pct,
                job_title = params$job_title
              )
            }
            
            invokeRestart("muffleMessage")
          }
        )
      )
      
      # Persist outputs
      writeLines(log_vec, file.path(jd, "log.txt"))
      saveRDS(as.data.frame(result_df), file.path(jd, "dice_result.rds"))
      saveRDS(modules$summary_df,       file.path(jd, "modules_summary.rds"))
      saveRDS(modules$membership_df,    file.path(jd, "modules_membership.rds"))
      saveRDS(modules$edges_by_module,  file.path(jd, "modules_edges.rds"))
      
      write_status(job_id, "finished", "Done.")
      list(job_id = job_id, df = result_df, log = log_vec, modules = modules)
      
    }, seed = TRUE) %...>% (function(res_list) {
      
      result_df      <- res_list$df
      dice_log_lines <- res_list$log
      modules_obj    <- res_list$modules
      
      if (!is.data.frame(result_df)) result_df <- as.data.frame(result_df)
      dice_result(result_df)
      
      if (length(dice_log_lines)) {
        add_log("----- DiCE console output -----")
        add_log(dice_log_lines)
      }
      phase_lines(grep("Genes in Phase", dice_log_lines, value = TRUE))
      
      modules_summary(modules_obj$summary_df)
      modules_membership(modules_obj$membership_df)
      modules_edges(modules_obj$edges_by_module)
      
      current_status("Finished.")
      removeModal()
      updateNavbarPage(session, "main_nav", "results")
      
    }) %...!% (function(e) {
      
      err_raw <- conditionMessage(e)
      fe <- friendly_dice_error(err_raw)
      
      write_status(
        job_id,
        state   = "error",
        message = clean_error_message(conditionMessage(e)),
        job_title = input$job_title
      )
      
      print(err_raw)
      
      add_log(paste("DiCE run failed:", conditionMessage(e)))
      current_status("Error.")
      removeModal()
      
      showModal(modalDialog(
        title = fe$title,
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$div(
          style = "font-size:15px; line-height:1.6;",
          
          ## --- Friendly explanation ---
          tags$p(tags$b("What went wrong:")),
          tags$p(fe$hint),
          
          tags$p(tags$b("Please verify the following before retrying:")),
          tags$ul(
            tags$li("Expression matrix: Genes should be provided as columns, samples as rows, and the class/phenotype label must be included in the last column. Sample names should not be placed in the first column."),
            tags$li("DGE file: must include gene name, logFC, p-value, and adjusted p-value."),
            tags$li("Treatment and control labels must exactly match values in the class column (case-sensitive).")
          ),
          
          tags$hr(),
          
          ## --- Technical details (collapsed) ---
          tags$details(
            tags$summary("Show technical details (for debugging)"),
            tags$pre(
              style = "white-space: pre-wrap; font-size: 12px; color: #444;",
              fe$details
            )
          )
        )
      ))
    })
  })
  
  
  
  ################################
  ## 15. Gene click -> boxplot
  ################################
  observeEvent(input$ppi_clicked_gene, {
    req(expr_df(), dge_df())
    
    g   <- input$ppi_clicked_gene
    dat <- expr_df()
    cond_col <- names(dat)[ncol(dat)]
    
    if (!(g %in% names(dat))) {
      showModal(modalDialog(
        title = paste("Gene:", g),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$em("This gene is not present as a column in the uploaded expression matrix.")
      ))
      return()
    }
    
    df_long <- data.frame(
      Condition  = dat[[cond_col]],
      Expression = dat[[g]]
    )
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Condition, y = Expression, fill = Condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
      ggplot2::labs(title = paste("Expression of", g), x = "", y = "Normalized expression") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "none")
    
    dge <- dge_df()
    gene_col <- dplyr::case_when(
      "Gene"      %in% names(dge) ~ "Gene",
      "gene"      %in% names(dge) ~ "gene",
      "Gene.Name" %in% names(dge) ~ "Gene.Name",
      TRUE ~ names(dge)[1]
    )
    row <- dge[dge[[gene_col]] == g, , drop = FALSE]
    
    if (nrow(row) == 0) {
      stats_html <- tags$em("No DGE statistics found for this gene.")
    } else {
      logFC <- if ("logFC"     %in% names(row)) signif(row$logFC[1], 4)     else NA
      pval  <- if ("P.Value"   %in% names(row)) signif(row$P.Value[1], 4)   else NA
      padj  <- if ("adj.P.Val" %in% names(row)) signif(row$adj.P.Val[1], 4) else NA
      
      stats_html <- tags$div(
        tags$h4("Differential Expression Statistics"),
        tags$p(HTML(paste0(
          "<b>logFC:</b> ", logFC, "<br>",
          "<b>P-value:</b> ", pval, "<br>",
          "<b>adj.P.Val:</b> ", padj
        )))
      )
    }
    
    showModal(modalDialog(
      title = paste("Gene:", g),
      size  = "l",
      easyClose = TRUE,
      footer = tagList(
        downloadButton("download_gene_plot", "Download PNG"),
        modalButton("Close")
      ),
      plotOutput("popup_gene_plot", height = "350px"),
      br(),
      stats_html
    ))
    
    output$popup_gene_plot <- renderPlot({ p })
  })
}

shinyApp(ui, server)
