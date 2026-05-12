about_tab <- function() {
  tabPanel(
    title = "About",
    value = "about",
    
    div(
      class = "page-container",
      
      # Hero section
      div(
        class = "hero-box",
        div(class = "hero-title", "About DiCE-Ex"),
        p(
          class = "hero-subtext",
          "DiCE-Ex is an interactive web platform for Differential Centrality-Ensemble (DiCE) analysis. ",
          "It helps identify biologically important genes by integrating differential expression, information gain, ",
          "condition-specific protein-protein interaction networks, and network centrality changes across biological states."
        )
      ),
      
      # Quick navigation
      div(
        class = "toc-box",
        h4(class = "subsection-title", "Quick navigation"),
        tags$ul(
          tags$li(tags$a(href = "#publication_section", "Publication")),
          tags$li(tags$a(href = "#docs_section", "Documentation")),
          tags$li(tags$a(href = "#workflow_section", "Workflow summary")),
          tags$li(tags$a(href = "#results_section", "Outputs and result interpretation")),
          tags$li(tags$a(href = "#modules_section", "DiCE module interpretation")),
          tags$li(tags$a(href = "#support_section", "Support"))
        )
      ),
      
      # Publication
      div(
        id = "publication_section",
        class = "info-card",
        h3(class = "section-title", "Publication"),
        p(tags$b("DiCE: differential centrality-ensemble analysis based on gene expression profiles and protein-protein interaction network")),
        p("Pashaei, E., Liu, S., Li, K.Y., Zang, Y., Yang, L., Lautenschlaeger, T., Huang, J., Lu, X., & Wan, J. (2025). DiCE: differential centrality-ensemble analysis based on gene expression profiles and protein–protein interaction network. Nucleic Acids Research, 53."),
        tags$a(
          "View article (Nucleic Acids Research)",
          href = "https://academic.oup.com/nar/article/53/13/gkaf609/8192812",
          target = "_blank",
          class = "btn-doc"
        )
      ),
      
      # GitHub
      div(
        class = "info-card",
        h3(class = "section-title", "Bug reports and feature requests"),
        p("Please use our GitHub issue tracker to report bugs, request new features, or ask questions about DiCE-Ex."),
        tags$a(
          "Report an Issue on GitHub",
          href = "https://github.com/Wanbioinfo/DiCE-Ex/issues",
          target = "_blank",
          class = "btn-doc"
        )
      ),
      
      # Documentation
      div(
        id = "docs_section",
        class = "info-card",
        h3(class = "section-title", "Documentation"),
        
        h4(class = "subsection-title", "Overview"),
        p(
          "DiCE-Ex implements the complete DiCE workflow, allowing users to upload bulk RNA-seq datasets and run all phases—from DGE-based candidate gene selection to IG filtering, network reconstruction, centrality analysis, and ensemble ranking."
        ),
        p(
          "The platform supports human and mouse datasets and provides interactive visualization of DiCE genes, PPI subnetworks, and community modules."
        ),
        
        br(),
        
        h4(class = "subsection-title", "Input requirements"),
        p("DiCE-Ex supports .csv, .tsv, .xlsx, and .rds files."),
        
        fluidRow(
          column(
            width = 6,
            div(
              class = "mini-card",
              h5("Normalized expression matrix"),
              p("Formatted as genes × samples. Do not place gene symbols as rownames. Place them under a column name `Gene`."),
            )
          ),
          column(
            width = 6,
            div(
              class = "mini-card",
              h5("Metadata file"),
              p("Formatted with samples as rows. Include sample IDs under a column named `sample_id` and phenotype/group labels under a column named `phenotype`.")
            )
          ),
          column(
            width = 6,
            div(
              class = "mini-card",
              h5("Differential expression results"),
              p("Expected columns: Gene.Symbol, P.Value, adj.P.Val, and logFC."),
              tags$img(src = "images/dge_image.png", class = "doc-img-sm")
            )
          )
        )
      ),
      
      # Workflow
      div(
        id = "workflow_section",
        class = "info-card",
        h3(class = "section-title", "Workflow Summary"),
        p("DiCE consists of five main phases that progressively refine and prioritize biologically meaningful genes."),
        
        div(
          class = "phase-box",
          div(class = "phase-title", "Phase I — Significance Filtering"),
          p("Constructs the initial candidate gene pool from user-supplied differential gene expression results using permissive statistical filtering."),
          tags$ul(
            tags$li(tags$b("P-value measure: "), "Adjusted p-value or raw P.Value."),
            tags$li(tags$b("Selected p-value threshold: "), "Default 0.05 to retain genes with potentially subtle but meaningful changes."),
            tags$li(tags$b("|log2FC| threshold: "), "Default 0 allows genes with small fold-changes to remain eligible.")
          )
        ),
        
        div(
          class = "phase-box",
          div(class = "phase-title", "Phase II — Information Gain (IG) Filtering"),
          p("Retains genes whose expression profiles are informative for distinguishing biological conditions."),
          tags$ul(
            tags$li(tags$b("IG method: "), "Information Gain (IG) , Weighted Information Gain (WIG) , or None. 
                    Standard IG measures how well each gene distinguishes the two conditions. 
                    WIG uses resampling to reduce class imbalance effects. 
                    Selecting None skips Phase II filtering and keeps all Phase I genes."),
            tags$li(tags$b("B for WIG: "), "Bootstrap iterations used to stabilize weighted information gain estimates."),
            tags$li(tags$b("IG cutoff rule: "), "Mean, median, IG > 0, or a custom threshold.")
          ),
        ),
        
        div(
          class = "phase-box",
          div(class = "phase-title", "Phase III — PPI Network Construction"),
          p("Builds condition-specific weighted PPI networks using STRING interactions and expression-derived correlations."),
          tags$ul(
            tags$li(tags$b("STRING database: "), "Known protein-protein interactions are restricted to candidate genes."),
            tags$li(tags$b("Condition-specific networks: "), "Separate networks are built for treatment and control."),
            tags$li(tags$b("Edge weights: "), "|correlation coefficient|; distance is defined as 1 - |correlation coefficient|."),
            tags$li(tags$b("Correlation type: "), "Pearson or Spearman.")
          )
        ),
        
        div(
          class = "phase-box",
          div(class = "phase-title", "Phase IV — Network Centrality and Ensemble Ranking"),
          p("Calculates multiple centrality metrics and quantifies shifts in network importance between biological conditions."),
          
          tags$ul(
            tags$li("Supported metrics include betweenness, eigenvector, closeness, PageRank, harmonic, strength, and authority."),
            tags$li("All genes are ranked based on the absolute differences in each selected centrality measure between the two conditions."),
            tags$li("A consensus ensemble rank is then computed using the ProductOfRank method.")
          ),
          
          div(
            class = "soft-note",
            tags$b("Note: "),
            "Select the most representative centrality measures, as highly correlated metrics may capture similar topological information. ",
            "For example, eigenvector centrality and authority often show strong correlations."
          ),
          
          fluidRow(
            column(
              6,
              tags$img(src = "images/tcga_prad_heatmap.png", class = "doc-img-sm"),
              p(class = "caption", "TCGA-PRAD dataset")
            ),
            column(
              6,
              tags$img(src = "images/nepc_heatmap.png", class = "doc-img-sm"),
              p(class = "caption", "NEPC 2019 dataset")
            )
          )
        ),
        
        div(
          class = "phase-box",
          div(class = "phase-title", "Phase V — DiCE Gene Identification"),
          p("Final DiCE genes are identified by evaluating the user-defined dice_rules and combining rule outcomes using the selected dice_logic."),
          tags$p(
            "A DiCE rule defines the criterion a gene must satisfy to be selected as a final DiCE gene. ",
            "For example, in ",
            tags$code("dice_centrality_rule(metric = 'betweenness', threshold_type = 'percent', threshold = 25)"),
            ", a gene passes if its betweenness centrality is within the top 25% in either the treatment or control network. ",
            "Multiple rules can be combined using ",
            tags$code("dice_logic"),
            ", where ",
            tags$code("'AND'"),
            " requires a gene to satisfy all rules and ",
            tags$code("'OR'"),
            " requires it to satisfy at least one rule."
          )
        )
      ),
      
      # Results
      div(
        id = "results_section",
        class = "info-card",
        h3(class = "section-title", "Outputs and Result Interpretation"),
        p("DiCE-Ex outputs are designed to support both gene-level prioritization and network-level interpretation."),
        
        h4(class = "subsection-title", "DiCE results table"),
        tags$table(
          class = "table table-bordered table-striped about-table",
          style = "max-width:1000px;",
          tags$thead(
            tags$tr(
              tags$th("Column"),
              tags$th("Description")
            )
          ),
          tags$tbody(
            tags$tr(tags$td("Gene.Symbol"), tags$td("Official gene symbol used as the primary identifier.")),
            tags$tr(tags$td("logFC"), tags$td("Log2 fold change between treatment and control groups.")),
            tags$tr(tags$td("P.Value"), tags$td("Raw p-value from differential expression analysis.")),
            tags$tr(tags$td("adj.P.Val"), tags$td("Adjusted p-value corrected for multiple testing.")),
            tags$tr(tags$td("IG / Weighted IG"), tags$td("Feature selection score showing how informative a gene is for separating classes.")),
            tags$tr(tags$td("Centrality measures"), tags$td("Network topology metrics computed in each condition-specific network.")),
            tags$tr(tags$td("Centrality differences"), tags$td("Absolute differences in centrality between treatment and control.")),
            tags$tr(tags$td("Centrality_rank columns"), tags$td("Ranks based on centrality differences; lower ranks indicate higher priority.")),
            tags$tr(tags$td("Pass_Rules"), tags$td("Lists the DiCE rules satisfied by the gene. For example, it may show which selected centrality-based/ensemble-based rules the gene passed.")),
            tags$tr(tags$td("Pass_Counts"), tags$td("Number of DiCE rules satisfied by the gene.")),
            tags$tr(tags$td("Phase"), tags$td("Workflow stage in which the gene was retained.")),
            tags$tr(tags$td("ProductOfRank"), tags$td("Combined score obtained by multiplying ranks across selected centralities.")),
            tags$tr(tags$td("Ensemble_Rank"), tags$td("Final consensus ranking; lower values indicate stronger candidates."))
          )
        ),
        
        tags$img(src = "images/dice_results.png", class = "doc-img"),
        p(class = "caption", "Example DiCE results table showing gene ranking and filtering options."),
        
        h4(class = "subsection-title", "How to interpret gene rankings"),
        tags$ul(
          tags$li("Top-ranked genes often combine expression relevance with strong changes in network influence."),
          tags$li("Some prioritized genes may show modest fold-change but large centrality shifts."),
          tags$li(tags$b("Phase = DiCE"), " indicates genes that passed the final centrality-based filtering step."),
          tags$li("Some Phase III genes may still rank highly because ranking and filtering capture different aspects of the analysis."),
          tags$li("Users are encouraged to inspect both top-ranked DiCE genes and highly ranked Phase III genes.")
        )
      ),
      
      # Modules
      div(
        id = "modules_section",
        class = "info-card",
        h3(class = "section-title", "Interpretation of DiCE Modules"),
        p("DiCE-Ex identifies community modules within the STRING-derived PPI network constructed from DiCE genes."),
        
        fluidRow(
          column(
            6,
            div(
              class = "mini-card",
              h5("Module summary"),
              p("Displays the detected modules, the number of nodes, and the number of edges within each module."),
              tags$img(src = "images/module_summary.png", class = "doc-img-sm"),
              p(class = "caption", "Example module summary table.")
            )
          ),
          column(
            6,
            div(
              class = "mini-card",
              h5("Module membership"),
              p("Lists genes assigned to each module along with within-module connectivity."),
              tags$img(src = "images/module_memb.png", class = "doc-img-sm"),
              p(class = "caption", "Example membership table.")
            )
          )
        ),
        
        fluidRow(
          column(
            6,
            div(
              class = "mini-card",
              h5("Network visualization"),
              p("Shows module-specific subnetworks with genes as nodes and PPIs as edges."),
              tags$img(src = "images/net_vis.png", class = "doc-img-sm"),
              p(class = "caption", "Example DiCE PPI module network.")
            )
          ),
          column(
            6,
            div(
              class = "mini-card",
              h5("Gene exploration"),
              p("Clicking a gene can display expression boxplots for treatment and control."),
              tags$img(src = "images/gene_explore.png", class = "doc-img-sm"),
              p(class = "caption", "Example gene exploration panel.")
            )
          )
        ),
        
        tags$ul(
          tags$li("Genes within the same module often participate in related functions or pathways."),
          tags$li("Users can perform GO or pathway enrichment analysis on module genes for biological interpretation.")
        )
      ),
      
      # Support
      div(
        id = "support_section",
        class = "info-card",
        h3(class = "section-title", "Support"),
        tags$ul(
          tags$li(
            tags$b("Report bugs / request features: "),
            tags$a(
              href = "https://github.com/Wanbioinfo/DiCE-Ex/issues",
              target = "_blank",
              "GitHub Issues"
            )
          ),
          tags$li(
            tags$b("Please include: "),
            "job title, dataset name (if public), a screenshot of the error, and the log output from the Results page."
          )
        )
      )
    )
  )
}
