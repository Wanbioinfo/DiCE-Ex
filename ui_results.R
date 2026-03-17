# ui_results.R

results_tab <- function() {
  tabPanel(
    "Results",
    value = "results",
    
    div(
      class = "page-container",
      
      br(),
      h3("DiCE Results"),
      
      uiOutput("job_status_ui"),
      uiOutput("error_help_ui"),
      uiOutput("job_progress_ui"),
      
      uiOutput("dice_job_title"),
      uiOutput("phase_summary"),
      br(),
      
      tabsetPanel(
        id   = "results_inner_tabs",
        type = "tabs",
        
        ## -------------------------------------------------
        ## 1. DiCE Results TAB
        ## -------------------------------------------------
        tabPanel(
          title = "DiCE Results",
          br(),
          
          # Top bar: download button (left) + phase + gene search (right)
          fluidRow(
            column(
              width = 6,
              uiOutput("download_buttons_ui")
            ),
            column(
              width = 6,
              div(
                style = "display:flex; justify-content:flex-end; align-items:center; gap:20px;",
                
                # Phase dropdown (drives input$dice_phase_filter)
                div(
                  style = "display:flex; align-items:center;",
                  tags$label("Phase:", style = "margin-right:8px; font-weight:600;"),
                  selectInput(
                    inputId  = "dice_phase_filter",
                    label    = NULL,
                    choices  = c("All phases" = "all"),
                    selected = "all",
                    width    = "150px"
                  )
                ),
                
                # Gene search box (drives input$dice_search)
                div(
                  style = "display:flex; align-items:center;",
                  tags$label("Gene:", style = "margin-right:8px; font-weight:600;"),
                  textInput(
                    inputId     = "dice_search",
                    label       = NULL,
                    placeholder = "e.g., CDK1",
                    width       = "180px"
                  )
                )
              )
            )
          ),
          
          br(),
          DT::DTOutput("dice_table")
        ),
        
        ## -------------------------------------------------
        ## 2. DiCE PPI modules TAB
        ## -------------------------------------------------
        tabPanel(
          "DiCE PPI modules",
          
          div(
            style = "margin-bottom:12px;",
            uiOutput("download_modules_ui")
          ),
          
          h3("Module statistics"),
          div(
            style = "max-width:720px;",
            DT::DTOutput("modules_stats_table")
          ),
          br(),
          
          h3("Module membership & module subnetwork"),
          helpText(
            "Use the module selector or the search box to filter rows. Click a gene row to visualize its subnetwork (within its module) on the right."
          ),
          helpText(
            "By clicking a node in the network, you can view the gene expression changes."
          ),
          
          fluidRow(
            column(
              width = 4,
              
              selectInput(
                "mm_module_filter",
                "Filter by module",
                choices  = c("All modules" = "all"),
                selected = "all",
                width = "70%"
              ),
              
              div(
                style = "margin:8px 0 10px 0;",
                textInput(
                  inputId     = "mm_gene_search",
                  label       = "Gene search:",
                  placeholder = "e.g., POLE2",
                  width       = "70%"
                )
              ),
              
              DT::DTOutput("modules_membership_table")
            ),
            
            column(
              width = 8,
              style = "padding-right:0;",
              
              div(
                class = "module-network-panel",
                visNetworkOutput(
                  "ppi_module_network",
                  height = "750px",
                  width  = "100%"
                )
              )
            )
          )
        )
      )
    )
  )
}
