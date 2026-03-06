home_tab <- function() {
  tabPanel(
    title = "Home",
    value = "home",
    
    fluidPage(
      div(
        class = "page-container",
        
        br(),
        
        ################################
        ## HERO SECTION
        ################################
        div(
          style = "
            border-radius:20px;
            background-color:#004b4b;
            padding:80px 40px;
            color:white;
            text-align:center;
            margin-bottom:40px;
          ",
          
          h1(
            "Welcome to DiCE-Ex",
            style = "color:white; font-weight:700; margin-bottom:15px;"
          ),
          
          p(
            "A Web Server for Network-Based Differential-Centrality Gene Analysis.",
            style = "font-size:1.2em; color:#e8f7f7; margin-bottom:30px;"
          ),
          
          actionButton(
            "go_run",
            label = "Run DiCE",
            class = "btn btn-light btn-lg",
            style = "padding:10px 30px; font-size:1.1em; border-radius:10px;"
          )
        ),
        
        ################################
        ## ABOUT SECTION
        ################################
        div(
          style = "margin-top:20px; margin-bottom:40px; text-align:center;",
          
          h2("About the Tool"),
          
          p(
            paste(
              "DiCE-Web provides an intuitive platform for differential gene discovery using the",
              "Differential Centrality-Ensemble (DiCE) framework.",
              "This multi-phase method integrates expression-based gene selection with",
              "network-topology analysis to identify disease-associated critical genes",
              "in an unbiased manner.",
              "Users can explore flexible options—including Information Gain (IG),",
              "multiple centrality metrics, and PPI subnetwork visualization—through",
              "a simple, interactive interface.",
              "The validated DiCE approach"
            ),
            tags$a(
              "(PMID: 40626556)",
              href = "https://pubmed.ncbi.nlm.nih.gov/40626556/",
              target = "_blank"
            ),
            paste(
              "has successfully identified key genes and pathways in prostate cancer",
              "and is applicable to both human and mouse transcriptomic datasets.",
              "DiCE-Web enables researchers to perform end-to-end gene discovery",
              "directly from their data with clarity and efficiency."
            ),
            style = "max-width:950px;
                    margin-left:auto;
                    margin-right:auto;
                    line-height:1.6;
                    color:#333;"
          )
        ),
        
        ################################
        ## GET STARTED SECTION
        ################################
        div(
          style = "margin-bottom:40px;",
          
          h2(
            HTML('<span style="color:#008b7c; font-weight:700;">Get Started</span> with DiCE-Ex'),
            align = "center"
          ),
          
          br(),
          
          fluidRow(
            column(
              width = 6,
              offset = 3,
              
              actionLink(
                inputId = "home_docs",
                label   = NULL,
                class   = "home-card-link",
                
                div(
                  class = "feature-card",
                  div(class = "feature-icon", shiny::icon("question-circle")),
                  h4("Documentation"),
                  p("See detailed guides, input formats, and workflow examples.")
                )
                
              )
            )
          )
        ),
        
        br(), br()
      )
    )
  )
}
