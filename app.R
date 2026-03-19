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

source("DiCE/calculate_NetCentralities.R")
source("DiCE/calculate_weightedIG.R")
source("DiCE/corr_calculations.R")
source("DiCE/ensemble_Ranking.R")
source("DiCE/createPPI.R")
source("DiCE/DiCE.R")
source("DiCE/fileReader.R")
source("DiCE/normalize_df_cols.R")
source("DiCE/Phase_1.R")
source("DiCE/Phase_2.R")
source("DiCE/Phase_3.R")
source("DiCE/Phase_4.R")
source("DiCE/Phase_5.R")
source("DiCE/protein_coding_filter.R")
source("DiCE/unweighted_module_analysis.R")
source("DiCE/weighted_module_analysis_helpers.R")
source("DiCE/weighted_module_analysis.R")

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

safe_input <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0) return(default)
  
  if (length(x) == 1) {
    if (is.na(x)) return(default)
    if (is.character(x) && !nzchar(trimws(x))) return(default)
  }
  
  x
}

jobs_root <- "jobs"
dir.create(jobs_root, showWarnings = FALSE, recursive = TRUE)

job_dir <- function(job_id) file.path(jobs_root, job_id)

print(job_dir)

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
      "Please ensure the differential gene expression analysis file includes a gene column (e.g., Gene.Name / Gene / SYMBOL / Gene.Symbol.),",
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
  if (is.null(job_id) || !nzchar(job_id)) return(NULL)
  
  f <- file.path(job_dir(job_id), "status.json")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  
  txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  if (length(txt) == 0) return(NULL)
  
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

      .page-container {
      max-width: 1200px;
      margin: 0 auto;
      padding: 25px 30px 60px 30px;
      }
  
      .hero-box {
        background: linear-gradient(135deg, #f7f9fc 0%, #eef4fb 100%);
        border: 1px solid #dbe6f3;
        border-radius: 18px;
        padding: 28px 30px;
        margin-bottom: 28px;
        box-shadow: 0 4px 14px rgba(0,0,0,0.05);
      }
  
      .hero-title {
        font-size: 32px;
        font-weight: 700;
        color: #1f3556;
        margin-bottom: 10px;
      }
  
      .hero-subtext {
        font-size: 17px;
        line-height: 1.7;
        color: #3b4c63;
        margin-bottom: 0;
      }
  
      .info-card {
        background: #ffffff;
        border: 1px solid #e4e9f0;
        border-radius: 16px;
        padding: 22px 24px;
        margin-bottom: 24px;
        box-shadow: 0 3px 10px rgba(0,0,0,0.04);
      }
  
      .section-title {
        font-size: 26px;
        font-weight: 700;
        color: #213a5b;
        margin-top: 0;
        margin-bottom: 12px;
      }
  
      .subsection-title {
        font-size: 21px;
        font-weight: 650;
        color: #29476d;
        margin-top: 16px;
        margin-bottom: 10px;
      }
  
      .phase-box {
        background: #fbfcfe;
        border-left: 5px solid #4c78a8;
        border-radius: 10px;
        padding: 16px 18px;
        margin-bottom: 20px;
      }
  
      .phase-title {
        font-size: 22px;
        font-weight: 700;
        color: #29476d;
        margin-bottom: 8px;
      }
  
      .soft-note {
        background: #f8fbff;
        border: 1px solid #d8e7f6;
        border-radius: 12px;
        padding: 14px 16px;
        margin: 12px 0 18px 0;
        color: #35506d;
      }
  
      .mini-card {
        background: #f9fbfd;
        border: 1px solid #e3ebf3;
        border-radius: 14px;
        padding: 16px;
        margin-bottom: 18px;
        height: 100%;
      }
  
      .mini-card h5 {
        margin-top: 0;
        font-weight: 700;
        color: #28476b;
      }
  
      .doc-img {
        width: 100%;
        max-width: 1000px;
        margin-top: 12px;
        border: 1px solid #d8dee8;
        border-radius: 10px;
        box-shadow: 0 3px 10px rgba(0,0,0,0.05);
      }
  
      .doc-img-sm {
        width: 100%;
        max-width: 520px;
        margin-top: 12px;
        border: 1px solid #d8dee8;
        border-radius: 10px;
        box-shadow: 0 3px 10px rgba(0,0,0,0.05);
      }
  
      .caption {
        font-size: 0.92em;
        color: #5a6675;
        margin-top: 8px;
      }
  
      .btn-doc {
        background-color: #2f6fb2;
        border-color: #2f6fb2;
        color: white !important;
        border-radius: 10px;
        padding: 10px 16px;
        font-weight: 600;
        text-decoration: none !important;
        display: inline-block;
        margin-top: 10px;
      }
  
      .btn-doc:hover {
        background-color: #24598f;
        border-color: #24598f;
        color: white !important;
      }
  
      .about-table th {
        background: #f3f7fb;
        color: #233b59;
      }
  
      .about-table td, .about-table th {
        vertical-align: top !important;
        padding: 10px 12px !important;
      }
    
      .toc-box {
        background: #ffffff;
        border: 1px solid #e4e9f0;
        border-radius: 14px;
        padding: 18px 20px;
        margin-bottom: 24px;
      }
  
      .toc-box ul {
        margin-bottom: 0;
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
  cleanup_jobs(keep_days = 5)
  
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
  
  output$phase5_selected_centralities_ui <- renderUI({
    cent_labels <- c(
      betweenness = "Betweenness",
      eigenvector = "Eigenvector",
      degree      = "Degree",
      pagerank    = "PageRank",
      closeness   = "Closeness",
      harmonic    = "Harmonic",
      authority   = "Authority",
      strength    = "Strength"
    )
    
    if (is.null(input$phase4_cents) || length(input$phase4_cents) == 0) {
      return(
        div(
          style = "
          border:1px solid #f0ad4e;
          background:#fffaf2;
          border-radius:8px;
          padding:10px 12px;
          margin-bottom:12px;
        ",
          "No centrality selected in Phase IV."
        )
      )
    }
    
    selected_labels <- unname(cent_labels[input$phase4_cents])
    
    div(
      style = "
      border:1px solid #d9edf7;
      background:#f7fcff;
      border-radius:8px;
      padding:10px 12px;
      margin-bottom:12px;
    ",
      tags$b("Selected centralities from Phase IV: "),
      paste(selected_labels, collapse = ", ")
    )
  })
  
  output$phase5_centrality_rules_ui <- renderUI({
    
    req(input$phase4_cents)
    
    cent_labels <- c(
      betweenness = "Betweenness",
      eigenvector = "Eigenvector",
      degree      = "Degree",
      pagerank    = "PageRank",
      closeness   = "Closeness",
      harmonic    = "Harmonic",
      authority   = "Authority",
      strength    = "Strength"
    )
    
    rule_ui_list <- lapply(input$phase4_cents, function(cent) {
      
      label <- cent_labels[[cent]]
      
      div(
        style = "
        border:1px solid #e9ecef;
        border-radius:8px;
        padding:12px;
        margin-bottom:12px;
        background:#fafafa;
      ",
        
        tags$h5(style = "margin-top:0; color:#004b4b;", label),
        
        fluidRow(
          column(
            6,
            selectInput(
              inputId = paste0("phase5_cutoff_type_", cent),
              label   = "Cutoff type",
              choices = c(
                "Mean"   = "mean",
                "Top K%" = "percent",
                "Top K"  = "rank"
              ),
              selected = "percent"
            )
          ),
          
          column(
            6,
            conditionalPanel(
              condition = sprintf("input['phase5_cutoff_type_%s'] == 'percent'", cent),
              numericInput(
                inputId = paste0("phase5_cutoff_percent_", cent),
                label   = "K (%) for top K% cutoff",
                value   = 25,
                min     = 1,
                max     = 100,
                step    = 1
              )
            ),
            
            conditionalPanel(
              condition = sprintf("input['phase5_cutoff_type_%s'] == 'rank'", cent),
              numericInput(
                inputId = paste0("phase5_cutoff_rank_", cent),
                label   = "K for top K cutoff",
                value   = 200,
                min     = 1,
                step    = 1
              )
            )
          )
        )
      )
    })
    
    do.call(tagList, rule_ui_list)
  })
  
  output$job_status_ui <- renderUI({
    jid <- job_id_qs()
    st  <- job_status()
    
    if (!is.null(jid) && is.null(st)) {
      return(
        div(
          class = "phase-summary",
          style = "border-left:5px solid #dc3545; background:#fff5f5;",
          tags$p(tags$b("This saved job is no longer available.")),
          tags$p("The job may have expired or been removed from server storage.")
        )
      )
    }
    
    if (is.null(st)) return(NULL)
    
    msg <- st$message
    show_msg <- is.character(msg) && length(msg) == 1 && !is.na(msg) && nzchar(msg) && !is_phase_message(msg)
    
    div(
      class = "phase-summary",
      tags$p(paste("Job:", st$job_id %||% "")),
      tags$p(paste("Status:", st$state %||% "unknown")),
      if (show_msg) tags$p(msg)
    )
  })
  
  output$dge_loaded_badge <- renderUI({
    req(dge_df())
    tags$div(
      style="margin-top:6px; color:#198754; font-size:13px;",
      icon("check-circle"), " Sample differential expression analysis data loaded"
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
  job_species_val <- reactiveVal(NULL)
  loaded_job_id   <- reactiveVal(NULL)
  
  modules_summary    <- reactiveVal(NULL)
  modules_membership <- reactiveVal(NULL)
  modules_edges      <- reactiveVal(NULL)
  modules_stats <- reactiveVal(NULL)
  modules_between_edges <- reactiveVal(NULL)
  all_module_edges <- reactiveVal(NULL)
  
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
        tags$b("Differential gene expression analysis file: "),
        "Must contain a valid gene identifier column (e.g., Gene or Gene.Name or Gene.Symbol.), ",
        "logFC, p-value, and adjusted p-value."
      ),
      tags$li(
        tags$b("Gene identifiers: "),
        "Must match between expression and Differential gene expression analysis files."
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
  
  observeEvent(list(job_status(), job_id_qs()), {
    st  <- job_status()
    jid <- job_id_qs()
    if (is.null(st) || is.null(jid)) return()
    
    updateNavbarPage(session, "main_nav", "results")
    
    if (!is.null(st$job_title) &&
        is.character(st$job_title) &&
        length(st$job_title) == 1 &&
        !is.na(st$job_title) &&
        nzchar(st$job_title)) {
      job_title_val(st$job_title)
    }
    
    if (identical(st$state, "finished") && !identical(loaded_job_id(), jid)) {
      jd <- job_dir(jid)
      
      f1 <- file.path(jd, "dice_result.rds")
      if (file.exists(f1)) dice_result(readRDS(f1))
      
      f2 <- file.path(jd, "modules_summary.rds")
      if (file.exists(f2)) modules_summary(readRDS(f2))
      
      f3 <- file.path(jd, "modules_membership.rds")
      if (file.exists(f3)) modules_membership(readRDS(f3))
      
      f4 <- file.path(jd, "modules_edges.rds")
      if (file.exists(f4)) modules_edges(readRDS(f4))
      
      f5 <- file.path(jd, "expr_input.rds")
      if (file.exists(f5)) expr_df(readRDS(f5))
      
      f6 <- file.path(jd, "dge_input.rds")
      if (file.exists(f6)) {
        dge_saved <- readRDS(f6)
        dge_saved <- as.data.frame(dge_saved)
        dge_saved <- normalize_dge_cols(dge_saved)
        dge_df(dge_saved)
      
        # overwrite old saved file so future loads are already fixed
        saveRDS(dge_saved, f6)
      }
      
      f7 <- file.path(jd, "module_stats.rds")
      if (file.exists(f7)) modules_stats(readRDS(f7))
      
      f8 <- file.path(jd, "between_module_edges.rds")
      if (file.exists(f8)) modules_between_edges(readRDS(f8))
      
      pf <- file.path(jd, "params.json")
      if (file.exists(pf)) {
        p <- jsonlite::read_json(pf, simplifyVector = TRUE)
        if (!is.null(p$species) && nzchar(p$species)) {
          job_species_val(tolower(p$species))
        }
      }
      
      f9 <- file.path(jd, "all_module_edges.rds")
      if (file.exists(f9)) all_module_edges(readRDS(f9))
      
      log_file <- file.path(jd, "log.txt")
      if (file.exists(log_file)) {
        log_lines <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) character(0))
        phase_lines(grep("Genes in Phase", log_lines, value = TRUE))
        dice_log(paste(log_lines, collapse = "\n"))
      }
      
      loaded_job_id(jid)
    }
  }, ignoreInit = FALSE)
  
  ################################
  ## 5. Upload paths (uploads OR sample files)
  ################################
  expr_path <- reactiveVal(NULL)
  dge_path  <- reactiveVal(NULL)
  
  observeEvent(input$expr_file, {
    req(input$expr_file)
    req(input$expr_file$datapath)
    req(nzchar(input$expr_file$datapath))
    
    path <- input$expr_file$datapath
    req(file.exists(path))
    
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
    req(input$dge_file$datapath)
    req(nzchar(input$dge_file$datapath))
    
    path <- input$dge_file$datapath
    req(file.exists(path))
    
    dge_path(path)
    
    ext <- tolower(tools::file_ext(input$dge_file$name))
    dge_tbl <- switch(
      ext,
      "csv"  = readr::read_csv(path, show_col_types = FALSE),
      "tsv"  = readr::read_tsv(path, show_col_types = FALSE),
      "txt"  = readr::read_tsv(path, show_col_types = FALSE),
      "xls"  = readxl::read_excel(path),
      "xlsx" = readxl::read_excel(path),
      "rds"  = {
        obj <- readRDS(path)
        if (is.matrix(obj)) obj <- as.data.frame(obj)
        if (!is.data.frame(obj)) stop("Differential gene expression analysis RDS must contain a data frame or matrix.")
        obj
      },
      readr::read_delim(path, delim = NULL, show_col_types = FALSE)
    )
    
    dge_tbl <- as.data.frame(dge_tbl)
    dge_tbl <- normalize_dge_cols(dge_tbl)
    dge_df(dge_tbl)
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
        "Gene.Symbol" %in% names(df) ~ "Gene.Symbol",
        "Gene.symbol" %in% names(df) ~ "Gene.symbol",
        "Gene.Name"   %in% names(df) ~ "Gene.Name",
        "Gene"        %in% names(df) ~ "Gene",
        "gene"        %in% names(df) ~ "gene",
        TRUE ~ NA_character_
      )
      if (!is.na(gene_col)) {
        keep <- grepl(q, as.character(df[[gene_col]]), ignore.case = TRUE)
        df   <- df[keep, , drop = FALSE]
      }
    }
    # Now add hyperlinks
    species_used <- job_species_val() %||% tolower(input$species %||% "human")
    df <- add_ncbi_gene_links(df, species = species_used)
    
    df
  })
  
  output$dice_table <- DT::renderDT({
    df <- dice_table_data()
    req(df)
  
    df_show <- df
  
    gene_col <- attr(df_show, "gene_col_linked")
    gene_idx <- if (!is.null(gene_col) && gene_col %in% names(df_show)) which(names(df_show) == gene_col) else NULL
  
    num_cols <- names(df_show)[vapply(df_show, is.numeric, logical(1))]
  
    # Adjust these if your column names differ
    sci_cols <- intersect(
      c("adj.P.Val", "P.Value", "adj.PVal", "PValue", "padj", "pvalue"),
      names(df_show)
    )
  
    int_cols <- intersect(
      c("Betweenness_rank", "EigenVector_rank", "Pass_Count", "Ensemble_Rank",
        "ProductRank", "pass_count", "Ensemble.Rank"),
      names(df_show)
    )
  
    fmt_sci <- function(x) {
      ifelse(is.na(x), NA_character_, format(x, scientific = TRUE, digits = 2, trim = TRUE))
    }
  
    fmt_int <- function(x) {
      ifelse(is.na(x), NA_character_, as.character(round(x)))
    }
  
    fmt_2dec <- function(x) {
      ifelse(is.na(x), NA_character_, sprintf("%.2f", x))
    }
  
    # First apply 2-decimal formatting to all numeric columns
    for (col in num_cols) {
      df_show[[col]] <- fmt_2dec(df_show[[col]])
    }
  
    # Override with scientific notation for p-value-like columns
    for (col in sci_cols) {
      df_show[[col]] <- fmt_sci(df[[col]])
    }
  
    # Override with integer formatting for rank/count columns
    for (col in int_cols) {
      df_show[[col]] <- fmt_int(df[[col]])
    }
  
    escape_cols <- if (!is.null(gene_idx)) setdiff(seq_along(df_show), gene_idx) else TRUE
  
    DT::datatable(
      df_show,
      filter = "none",
      rownames = FALSE,
      escape = escape_cols,
      options = list(
        dom = "lrtip",
        pageLength = 20,
        scrollX = TRUE
      ),
      selection = "none"
    )
  })
  
  output$modules_stats_table <- DT::renderDT({
    df <- modules_stats()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    
    df <- df %>%
      dplyr::rename(
        `Number of Nodes` = Num_Nodes,
        `Number of Edges` = Num_Edges
      )
    
    DT::datatable(
      df,
      filter = "none",
      rownames = FALSE,
      options = list(
        dom = "tip",
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = FALSE,
        columnDefs = list(
          list(width = "120px", targets = 0),
          list(width = "120px", targets = 1),
          list(width = "120px", targets = 2)
        )
      )
    )
  })
  
  output$modules_between_edges_table <- DT::renderDT({
    df <- modules_between_edges()
    req(df)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    
    df <- df %>%
      dplyr::rename(
        A = Module_A,
        B = Module_B,
        `edge_count (inter M)` = Num_Edges
      )
    
    DT::datatable(
      df,
      filter = "none",
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10, dom = "lrtip")
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
      df <- df[df$Module == input$mm_module_filter, , drop = FALSE]
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
      s  <- modules_summary()
      ms <- modules_stats()
      be <- modules_between_edges()
      m  <- modules_membership()
      ae <- all_module_edges()
      req(s, ms, be, m, ae)
      
      wb <- openxlsx::createWorkbook()
      
      openxlsx::addWorksheet(wb, "Module summary")
      openxlsx::writeData(wb, "Module summary", as.data.frame(s))
      
      openxlsx::addWorksheet(wb, "Module statistics")
      openxlsx::writeData(wb, "Module statistics", as.data.frame(ms))
      
      openxlsx::addWorksheet(wb, "Between-module edges")
      openxlsx::writeData(wb, "Between-module edges", as.data.frame(be))
      
      openxlsx::addWorksheet(wb, "Module membership")
      openxlsx::writeData(wb, "Module membership", as.data.frame(m))
      
      openxlsx::addWorksheet(wb, "All edges")
      openxlsx::writeData(wb, "All edges", as.data.frame(ae))
      
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
    input$mm_module_filter
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
  
  output$module_overview_network <- renderVisNetwork({
    ms <- modules_stats()
    be <- modules_between_edges()
    req(ms)
    
    if (!is.data.frame(ms)) ms <- as.data.frame(ms)
    if (!is.null(be) && !is.data.frame(be)) be <- as.data.frame(be)
    
    if (nrow(ms) == 0) {
      return(visNetwork(data.frame(), data.frame()) %>% visOptions(nodesIdSelection = FALSE))
    }
    
    mod_col   <- if ("Module" %in% names(ms)) "Module" else names(ms)[1]
    nodes_col <- if ("Num_Nodes" %in% names(ms)) "Num_Nodes" else names(ms)[2]
    intra_col <- if ("Num_Edges" %in% names(ms)) "Num_Edges" else names(ms)[3]
    
    nodes_df <- data.frame(
      id    = as.character(ms[[mod_col]]),
      label = as.character(ms[[mod_col]]),
      value = as.numeric(ms[[nodes_col]]),
      title = paste0(
        "<b>", ms[[mod_col]], "</b><br>",
        "Number of nodes: ", ms[[nodes_col]], "<br>",
        "Number of edges: ", ms[[intra_col]]
      ),
      stringsAsFactors = FALSE
    )
    all_mods <- sort(unique(as.character(modules_membership()$Module)))
    mod_cols <- grDevices::hcl.colors(length(all_mods), palette = "Set 2")
    names(mod_cols) <- all_mods
    
    if (!is.null(be) && nrow(be) > 0) {
      modA_col <- if ("Module_A" %in% names(be)) "Module_A" else names(be)[1]
      modB_col <- if ("Module_B" %in% names(be)) "Module_B" else names(be)[2]
      be_col   <- if ("Num_Edges" %in% names(be)) "Num_Edges" else names(be)[3]
      
      edges_df <- data.frame(
        from  = as.character(be[[modA_col]]),
        to    = as.character(be[[modB_col]]),
        value = as.numeric(be[[be_col]]),
        title = paste0("Between-module edges: ", be[[be_col]]),
        stringsAsFactors = FALSE
      )
    } else {
      edges_df <- data.frame(from = character(0), to = character(0), value = numeric(0))
    }
    
    mod_levels <- nodes_df$id
    pal <- grDevices::hcl.colors(length(mod_levels), "Set 2")
    names(pal) <- mod_levels
    
    nodes_df$color.background <- mod_cols[nodes_df$Module]
    nodes_df$color.border     <- mod_cols[nodes_df$Module]
    
    visNetwork(nodes_df, edges_df, width = "100%", height = "430px") %>%
      visNodes(
        shape = "dot",
        scaling = list(min = 20, max = 70),
        font = list(size = 24, face = "bold")
      ) %>%
      visEdges(
        smooth = FALSE,
        scaling = list(min = 1, max = 12)
      ) %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, hover = TRUE),
        nodesIdSelection = FALSE
      ) %>%
      visIgraphLayout(layout = "layout_with_fr") %>%
      visPhysics(enabled = FALSE)
  })
  
  output$ppi_module_network <- renderVisNetwork({
    nodes <- module_nodes()
    edges <- module_edges()
    req(nodes, edges)
    
    nodes_df <- data.frame(
      id     = nodes$Gene,
      label  = nodes$Gene,
      title  = paste0(
        nodes$Gene,
        "<br>Module: ", nodes$Module,
        "<br>Degree: ", nodes$Degree_inModule
      ),
      value  = nodes$Degree_inModule,
      Module = as.character(nodes$Module),
      stringsAsFactors = FALSE
    )
    
    all_mods <- sort(unique(as.character(modules_membership()$Module)))
    mod_cols <- grDevices::hcl.colors(length(all_mods), palette = "Set 2")
    names(mod_cols) <- all_mods
    
    edges_df <- data.frame(
      from = edges$Gene1,
      to   = edges$Gene2,
      stringsAsFactors = FALSE
    )
    
    nodes_df <- nodes_df[!is.na(nodes_df$id) & nzchar(nodes_df$id), , drop = FALSE]
    edges_df <- edges_df[
      !is.na(edges_df$from) & nzchar(edges_df$from) &
        !is.na(edges_df$to) & nzchar(edges_df$to),
      , drop = FALSE
    ]
    
    if (nrow(nodes_df) <= 1 || nrow(edges_df) == 0) {
      return(
        visNetwork(data.frame(), data.frame()) %>%
          visOptions(nodesIdSelection = FALSE)
      )
    }
    
    nodes_df$color.background <- mod_cols[nodes_df$Module]
    nodes_df$color.border     <- mod_cols[nodes_df$Module]
    
    highlight_gene <- selected_gene_id()
    
    if (!is.null(highlight_gene) && nzchar(highlight_gene)) {
      keep_edges <- edges_df$from == highlight_gene | edges_df$to == highlight_gene
      edges_sub  <- edges_df[keep_edges, , drop = FALSE]
      
      if (nrow(edges_sub) > 0) {
        keep_nodes <- unique(c(edges_sub$from, edges_sub$to))
        nodes_df <- nodes_df[nodes_df$id %in% keep_nodes, , drop = FALSE]
        edges_df <- edges_sub
      }
      
      nodes_df$color.background <- ifelse(
        nodes_df$id == highlight_gene,
        "#e74c3c",
        nodes_df$color.background
      )
      nodes_df$color.border <- ifelse(
        nodes_df$id == highlight_gene,
        "#c0392b",
        nodes_df$color.border
      )
    }
    
    max_nodes <- 200
    if (is.null(highlight_gene) && nrow(nodes_df) > max_nodes) {
      top_ids  <- nodes_df$id[order(-nodes_df$value)][1:max_nodes]
      nodes_df <- nodes_df[nodes_df$id %in% top_ids, , drop = FALSE]
      edges_df <- edges_df[
        edges_df$from %in% nodes_df$id & edges_df$to %in% nodes_df$id,
        , drop = FALSE
      ]
    }
    
    file_name <- if (!is.null(highlight_gene) && nzchar(highlight_gene)) {
      paste0(highlight_gene, "_module_network")
    } else {
      paste0("Module_", selected_module_id(), "_network")
    }
    
    visNetwork(nodes_df, edges_df, width = "100%", height = "780px") %>%
      visNodes(
        shape = "dot",
        size = 32,
        font = list(size = 28)
      ) %>%
      visEdges(
        smooth = FALSE,
        color = list(opacity = 0.75)
      ) %>%
      visOptions(
        highlightNearest = FALSE,
        nodesIdSelection = FALSE
      ) %>%
      visIgraphLayout(layout = "layout_with_fr", 
                      randomSeed = 123,
                      type = "full") %>%
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
      
      updateCheckboxGroupInput(
        session,
        "phase4_cents",
        selected = c("betweenness", "eigenvector")
      )
      
      sample_expr <- "sample_data/sample_Human_data_geneExp.csv"
      sample_dge  <- "sample_data/sample_Human_data_DGE.csv"
      
      req(file.exists(sample_expr), file.exists(sample_dge))
      
      expr_path(sample_expr)
      dge_path(sample_dge)
      
      expr_df(vroom::vroom(sample_expr, show_col_types = FALSE))
      dge_tbl <- vroom::vroom(sample_dge, show_col_types = FALSE)
      dge_tbl <- as.data.frame(dge_tbl)
      dge_tbl <- normalize_dge_cols(dge_tbl)
      dge_df(dge_tbl)
      
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
          "Human sample expression and Differential expression analysis files have been loaded successfully.<br><br>",
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
      
      updateCheckboxGroupInput(
        session,
        "phase4_cents",
        selected = c("betweenness", "eigenvector")
      )
      
      sample_expr <- "sample_data/sample_Mouse_data_geneExp.csv"
      sample_dge  <- "sample_data/sample_Mouse_data_DGE.csv"
      
      expr_path(sample_expr)
      dge_path(sample_dge)
      
      expr_df(vroom::vroom(sample_expr, show_col_types = FALSE))
      
      dge_tbl <- vroom::vroom(sample_dge, show_col_types = FALSE)
      dge_tbl <- as.data.frame(dge_tbl)
      dge_tbl <- normalize_dge_cols(dge_tbl)
      dge_df(dge_tbl)
      
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
          "Mouse sample expression and differential expression analysis files have been loaded successfully.<br><br>",
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
    
    req(input$phase4_cents)
    shiny::validate(
      shiny::need(length(input$phase4_cents) > 0, "Please select at least one centrality in Phase IV.")
    )
    
    dice_rules <- list()
    
    selection_mode <- input$phase5_selection_mode %||% "both"
    
    if (selection_mode %in% c("centrality_only", "both")) {
      
      req(input$phase4_cents)
      shiny::validate(
        shiny::need(length(input$phase4_cents) > 0,
                    "Please select at least one centrality in Phase IV.")
      )
      
      centrality_rules <- lapply(input$phase4_cents, function(cent) {
        
        cutoff_type <- safe_input(input[[paste0("phase5_cutoff_type_", cent)]], "percent")
        
        if (identical(cutoff_type, "mean")) {
          dice_centrality_rule(
            metric = cent,
            threshold_type = "mean"
          )
          
        } else if (identical(cutoff_type, "percent")) {
          dice_centrality_rule(
            metric = cent,
            threshold_type = "percent",
            threshold = safe_input(input[[paste0("phase5_cutoff_percent_", cent)]], 25)
          )
          
        } else if (identical(cutoff_type, "rank")) {
          dice_centrality_rule(
            metric = cent,
            threshold_type = "rank",
            threshold = safe_input(input[[paste0("phase5_cutoff_rank_", cent)]], 200)
          )
          
        } else {
          stop("Invalid Phase V cutoff type for centrality: ", cent)
        }
      })
      
      dice_rules <- c(dice_rules, centrality_rules)
    }
    
    if (selection_mode %in% c("ensemble_only", "both")) {
      
      ensemble_rank_cutoff <- input$phase5_ensemble_rank
      if (is.null(ensemble_rank_cutoff) || length(ensemble_rank_cutoff) == 0 || is.na(ensemble_rank_cutoff)) {
        ensemble_rank_cutoff <- 200
      }
      
      dice_rules <- c(
        dice_rules,
        list(
          dice_ensemble_rule(
            threshold_type = "rank",
            threshold = ensemble_rank_cutoff
          )
        )
      )
    }
    
    if (isTRUE(input$phase5_use_ensemble_rule)) {
      
      ensemble_rank_cutoff <- input$phase5_ensemble_rank
      if (is.null(ensemble_rank_cutoff) || length(ensemble_rank_cutoff) == 0 || is.na(ensemble_rank_cutoff)) {
        ensemble_rank_cutoff <- 200
      }
      
      dice_rules <- c(
        dice_rules,
        list(
          dice_ensemble_rule(
            threshold_type = "rank",
            threshold = ensemble_rank_cutoff
          )
        )
      )
    }
    
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
      dice_rules   = dice_rules,
      dice_logic   = input$phase5_dice_logic
    )
    
    jsonlite::write_json(
      params,
      path = file.path(jd, "params.json"),
      auto_unbox = TRUE,
      pretty = TRUE
    )
    
    dge_file  <- dge_path()
    expr_file <- expr_path()
    
    req(expr_file, dge_file)
    req(nzchar(expr_file), nzchar(dge_file))
    req(file.exists(expr_file), file.exists(dge_file))
    
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
              corr_pval_cutoff      = 1,
              centrality_list       = params$centralities,
              dice_rules            = params$dice_rules,
              dice_logic            = params$dice_logic,
              remove_pc_genes  = FALSE
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
      saveRDS(modules$module_stats_df,          file.path(jd, "module_stats.rds"))
      saveRDS(modules$between_module_edges_df,  file.path(jd, "between_module_edges.rds"))
      saveRDS(modules$all_edges_df, file.path(jd, "all_module_edges.rds"))
      
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
      modules_stats(modules_obj$module_stats_df)
      modules_between_edges(modules_obj$between_module_edges_df)
      all_module_edges(modules_obj$all_edges_df)
      
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
            tags$li("Differential gene expression analysis file: must include gene name, logFC, p-value, and adjusted p-value."),
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
    g <- input$ppi_clicked_gene
    if (is.null(g) || !nzchar(g)) return()
    
    dat <- expr_df()
    dge <- dge_df()
    
    if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0 || ncol(dat) < 2) {
      showModal(modalDialog(
        title = paste("Gene:", g),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$em("Expression data are not available for this bookmarked job.")
      ))
      return()
    }
    
    if (is.null(dge) || !is.data.frame(dge) || nrow(dge) == 0) {
      showModal(modalDialog(
        title = paste("Gene:", g),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$em("Differential expression data are not available for this bookmarked job.")
      ))
      return()
    }
    
    cond_col <- names(dat)[ncol(dat)]
    if (is.null(cond_col) || !nzchar(cond_col)) {
      showModal(modalDialog(
        title = paste("Gene:", g),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$em("Could not identify the condition/class column in the expression matrix.")
      ))
      return()
    }
    
    if (!(g %in% names(dat))) {
      showModal(modalDialog(
        title = paste("Gene:", g),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$em("This gene is not present as a column in the saved expression matrix.")
      ))
      return()
    }
    
    df_long <- data.frame(
      Condition = dat[[cond_col]],
      Expression = dat[[g]],
      stringsAsFactors = FALSE
    )
    
    df_long <- df_long[!is.na(df_long$Expression), , drop = FALSE]
    
    if (nrow(df_long) == 0) {
      showModal(modalDialog(
        title = paste("Gene:", g),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$em("No expression values available for this gene.")
      ))
      return()
    }
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Condition, y = Expression, fill = Condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
      ggplot2::labs(title = paste("Expression of", g), x = "", y = "Normalized expression") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(legend.position = "none")
    
    gene_col <- dplyr::case_when(
      "Gene.Symbol" %in% names(dge) ~ "Gene.Symbol",
      "Gene.symbol" %in% names(dge) ~ "Gene.symbol",
      "Gene.Name"   %in% names(dge) ~ "Gene.Name",
      "Gene"        %in% names(dge) ~ "Gene",
      "gene"        %in% names(dge) ~ "gene",
      TRUE ~ NA_character_
    )
    
    if (is.na(gene_col)) {
      stats_html <- tags$em("No recognizable gene column found in the saved differential gene expression analysis table.")
    } else {
      row <- dge[as.character(dge[[gene_col]]) == g, , drop = FALSE]
      
      if (nrow(row) == 0) {
        stats_html <- tags$em("No differential gene expression analysis statistics found for this gene.")
      } else {
        logFC <- if ("logFC" %in% names(row)) {
          sprintf("%.2f", row$logFC[1])
        } else {
          NA
        }
        
        pval <- if ("P.Value" %in% names(row)) {
          format(row$P.Value[1], scientific = TRUE, digits = 2, trim = TRUE)
        } else if ("pvalue" %in% names(row)) {
          format(row$pvalue[1], scientific = TRUE, digits = 2, trim = TRUE)
        } else {
          NA
        }
        
        padj <- if ("adj.P.Val" %in% names(row)) {
          format(row$adj.P.Val[1], scientific = TRUE, digits = 2, trim = TRUE)
        } else if ("adjp" %in% names(row)) {
          format(row$adjp[1], scientific = TRUE, digits = 2, trim = TRUE)
        } else {
          NA
        }
        
        stats_html <- tags$div(
          tags$h4("Differential Expression Statistics"),
          tags$p(HTML(paste0(
            "<b>logFC:</b> ", logFC, "<br>",
            "<b>P-value:</b> ", pval, "<br>",
            "<b>adj.P.Val:</b> ", padj
          )))
        )
      }
    }
    
    output$popup_gene_plot <- renderPlot({ p })
    
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
  })
}

shinyApp(ui, server)
