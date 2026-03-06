# ui_run.R

run_tab <- function() {
  tabPanel(
    "Run DiCE",
    value = "run",
    fluidPage(
      div(
        class = "page-container",   # keeps margins consistent with home page
        
        br(),
        
        h2("Run DiCE"),
        br(),
        
        fluidRow(
          column(
            width = 4,
            textInput(
              "job_title",
              label = "Job title"
            )
          ),
          
          column(
            width = 8,
            br(),
            
            # ---- Align buttons RIGHT ----
            div(
              style = "display:flex; justify-content:flex-end; gap:12px; align-items:center;",
              
              div(
                class = "sample-select-wrap",
                selectInput(
                  "sample_dataset",
                  label = NULL,
                  choices = c(
                    "Select sample dataset" = "",
                    "Human sample data" = "human",
                    "Mouse sample data" = "mouse"
                  ),
                  selected = "",
                  width = "220px",
                  selectize = FALSE
                )
              ),
              
              downloadButton(
                "download_sample_zip",
                "Download sample data (.zip)",
                class = "btn btn-success"
              )
            )
          )
        ),
        
        h3("Input data & parameters"),
        
        # ---- Data uploads ----
        fluidRow(
          
          column(
            width = 6,
            column(
              width = 8,
              fileInput(
                "expr_file",
                label = "Gene expression data",
                accept = c(".csv", ".tsv", ".txt", ".xlsx",".rds")
              ),
              uiOutput("expr_loaded_badge"),
            ),
            column(
              width = 12,
              helpText(
                "Upload a normalized expression matrix with genes as columns, ",
                "samples as rows, and the phenotype/class label in the last column."
              )
            )
          ),
          
          column(
            width = 6,
            column(
              width = 8,
              fileInput(
                "dge_file",
                label = "Differential gene expression results",
                accept = c(".csv", ".tsv", ".txt", ".xlsx",".rds")
              ),
              uiOutput("dge_loaded_badge"),
            ),
            column(
              width = 12,
              helpText(
                "Upload a table with gene identifiers, logFC, p-value, ",
                "and adjusted p-value."
              )
            )
          )
        ),
        
        tags$hr(),
        
        # ---- Group labels: treatment vs control ----
        fluidRow(
          column(
            width = 6,
            textInput(
              "group_treat",
              label = "Treatment / case group label",
              placeholder = "e.g., AD, Tumor, Treatment"
            ),
            helpText(
              "This must match one of the values in the phenotype/class column ",
              "of your expression matrix."
            )
          )
        ),
        fluidRow(
          column(
            width = 6,
            textInput(
              "group_control",
              label = "Control / reference group label",
              placeholder = "e.g., ND, Normal, Baseline"
            ),
            helpText(
              "This must match the control group label in the phenotype/class column."
            )
          )
        ),
        
        tags$hr(),
        
        # ---- Species ----
        fluidRow(
          column(
            width = 12,
            h4("Species"),
            radioButtons(
              "species",
              label = NULL,
              choices = c("Human", "Mouse"),
              inline  = TRUE,
              selected = "Human"
            )
          )
        ),
        
        tags$hr(),
        
        # ---- Parameters by phase ----
        h4("Parameters by phase"),
        
        # Phase I – Filter genes based on Significance metrics 
        
        div(
          style = "
            border:1px solid #e5e5e5;
            border-radius:10px;
            margin-bottom:12px;
            padding:0;
          ",
          
          # --- Toggle Header with Arrow ---
          tags$div(
            class = "toggle-header",
            `data-toggle` = "collapse",
            `data-target` = "#phase1_panel",
            style = "
              padding:12px 15px;
              background-color:#f5f5f5;
              border-radius:10px 10px 0 0;
              cursor:pointer;
              font-weight:600;
              color:#004b4b;
              position:relative;
            ",
            "Phase I – Significance Filtering",
            tags$span(class="toggle-arrow")
          ),
          
          # --- Collapsible Content (closed by default) ---
          div(
            id = "phase1_panel",
            class = "collapse",   # <--- collapsed by default
            style = "padding:15px;",
            
            p(
              style = "margin-bottom:15px;",
              "Filter genes based on differential expression analysis (DEA) results."
            ),
            
            fluidRow(
              column(
                4,
                selectInput(
                  "significant_metric",
                  "P-value measure",
                  choices = c(
                    "Adjusted p-value" = "adj.P.Val",
                    "P-value"          = "P.Value"
                  ),
                  selected = "adj.P.Val"
                )
              ),
              column(
                4,
                numericInput(
                  "significant_thresh",
                  "Threshold for selected p-value",
                  value = 0.05,
                  min   = 0,
                  max   = 1,
                  step  = 0.001
                )
              ),
              column(
                4,
                numericInput(
                  "phase1_logfc_thresh",
                  "Threshold for |log2FC|",
                  value = 0,
                  step  = 0.1
                )
              )
            )
          )
        ) ,
        
        # Phase II – Feature filtering
        div(
          style = "
            border:1px solid #e5e5e5;
            border-radius:10px;
            margin-bottom:12px;
            padding:0;
          ",
          
          # ---- Toggle Header ----
          tags$div(
            class = "toggle-header collapsed",
            `data-toggle` = "collapse",
            `data-target` = "#phase2_panel",
            role = "button",
            `aria-expanded` = "false",
            style = "
              padding:12px 15px;
              background-color:#f5f5f5;
              border-radius:10px 10px 0 0;
              cursor:pointer;
              font-weight:600;
              color:#004b4b;
              position:relative;
              margin-bottom:0;
            ",
            "Phase II – IG-based Feature Selection",
            tags$span(class="toggle-arrow")
          ),
          
          # ---- Collapsible Content ----
          div(
            id = "phase2_panel",
            class = "collapse",
            style = "padding:15px;",
            
            p(
              style = "margin-bottom:10px;",
              "Identify discriminative genes using Information Gain (IG) or Weighted IG (WIG) or you can skip this phase."
            ),
            
            fluidRow(
              column(
                4,
                selectInput(
                  "phase2_ig_method",
                  "Information Gain method",
                  choices = c(
                    "Information Gain (IG)"           = "IG",
                    "Weighted Information Gain (WIG)" = "wIG",
                    "None" = "none"
                  ),
                  selected = "IG"
                )
              ),
              
              column(
                4,
                numericInput(
                  "phase2_B",
                  "B (bootstrap resamples for WIG)",
                  value = 300,
                  min   = 10,
                  step  = 10
                )
              ),
              
              column(
                4,
                selectInput(
                  "phase2_ig_cutoff",
                  "IG cutoff rule",
                  choices = c(
                    "IG > mean (all genes)" = "all_mean",
                    "IG > median (all genes)" = "all_median",
                    "IG > mean (IG > 0)" = "nonzero_mean",
                    "IG > median (IG > 0)" = "nonzero_median",
                    "IG > 0" = "all_nonzero",
                    "Custom" = "custom"
                  ),
                  selected = "all_mean"
                ),
                
                conditionalPanel(
                  condition = "input.phase2_ig_cutoff == 'custom'",
                  numericInput(
                    "phase2_custom_ig",
                    "Enter custom IG cutoff",
                    value = 0.01,
                    min = 0,
                    step = 0.001
                  )
                )
              )
            )
          )
        ),
        
        # Phase III – Network construction
        div(
          style = "
            border:1px solid #e5e5e5;
            border-radius:10px;
            margin-bottom:12px;
            padding:0;
          ",
          
          tags$div(
            class = "toggle-header collapsed",
            `data-toggle` = "collapse",
            `data-target` = "#phase3_panel",
            role = "button",
            `aria-expanded` = "false",
            style = "
              padding:12px 15px;
              background-color:#f5f5f5;
              border-radius:10px 10px 0 0;
              cursor:pointer;
              font-weight:600;
              color:#004b4b;
              position:relative;
              margin-bottom:0;
            ",
            "Phase III – PPI Network Construction",
            tags$span(class = "toggle-arrow")
          ),
          
          # ---- Collapsible Body ----
          div(
            id = "phase3_panel",
            class = "collapse",
            style = "padding:15px;",
            
            p(
              style = "margin-bottom:10px;",
              "Build condition-specific STRING PPI networks with correlation-based edge weights."
            ),
            
            fluidRow(
              column(
                4,
                selectInput(
                  "phase3_corr",
                  "Correlation type",
                  choices = c("Pearson", "Spearman"),
                  selected = "Pearson"
                )
              )
            )
          )
        ),
        
        # Phase IV – Centrality aggregation
        div(
          style = "
            border:1px solid #e5e5e5;
            border-radius:10px;
            margin-bottom:12px;
            padding:0;
          ",
          
          # ---- Toggle Header ----
          tags$div(
            class = "toggle-header collapsed",
            `data-toggle` = "collapse",
            `data-target` = "#phase4_panel",
            role = "button",
            `aria-expanded` = "false",
            style = "
              padding:12px 15px;
              background-color:#f5f5f5;
              border-radius:10px 10px 0 0;
              cursor:pointer;
              font-weight:600;
              color:#004b4b;
              position:relative;
              margin-bottom:0;
            ",
            "Phase IV – Centrality-based DiCE Filtering",
            tags$span(class = "toggle-arrow")
          ),
          
          # ---- Collapsible Body ----
          div(
            id = "phase4_panel",
            class = "collapse",
            style = "padding:15px;",
                    
            p(
              style = "margin-bottom:8px;",
              "Perform network topology analysis and filter DiCE genes using centrality metrics."
            ),
            
            fluidRow(
              # Centrality measures
              column(
                6,
                checkboxGroupInput(
                  "phase4_cents",
                  "Centrality measures",
                  choices = c(
                    "Betweenness"   = "betweenness",
                    "Eigenvector"   = "eigenvector",
                    "Degree"        = "degree",
                    "PageRank"      = "pagerank",
                    "Closeness"     = "closeness",
                    "Harmonic"      = "harmonic",
                    "Authority"     = "authority",
                    "Strength"      = "strength"
                  ),
                  selected = c("betweenness", "eigenvector")
                )
              ),
              
              # DiCE cutoff rule + K%
              column(
                6,
                selectInput(
                  "phase4_dice_cutoff",
                  "DiCE cutoff",
                  choices = c(
                    "Mean"    = "mean",
                    "Median"   = "median",
                    "Top K%"   = "topK"
                  ),
                  selected = "topK"
                ),
                conditionalPanel(
                  "input.phase4_dice_cutoff == 'topK'",
                  numericInput(
                    "phase4_topK",
                    "K (%) for top K% cutoff",
                    value = 25,
                    min   = 1,
                    max   = 100,
                    step  = 1
                  )
                )
              )
            )
          )
        ),
        
        tags$hr(),
        
        div(
          style = "text-align:right; padding-bottom:60px;",
          actionButton(
            "run_dice",
            "Run DiCE",
            class = "btn btn-success"
          )
        )
      )
    )
  )
}
