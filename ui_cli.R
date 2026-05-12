# ui_cli.R
cli_tab <- function() {
  tabPanel(
    title = "DiCE CLI",
    value = "cli",
    
    div(
      class = "page-container",
      
      # Header
      tags$div(
        style = "margin-top: 10px; margin-bottom: 18px;",
        tags$h2("DiCE CLI", style = "margin-bottom:6px;"),
        tags$p(
          style = "color:#555; font-size:15px; line-height:1.6;",
          "Run DiCE from the command line using a JSON config file. ",
          "This is useful for reproducible runs, scripting, and running on servers/HPC."
        )
      ),
      
      tags$hr(),
      
      tags$h3("Source Code"),
      
      tags$p(
        "The DiCE source code and CLI scripts are available on GitHub."
      ),
      
      tags$a(
        href = "https://github.com/Wanbioinfo/DiCE-Ex",
        target = "_blank",
        class = "btn btn-success",
        icon("github"),
        "Access Source Code on GitHub"
      ),
      
      tags$h3("Environment Setup"),
      
      tags$p(
        "Before running the DiCE CLI, install the required R packages."
      ),
      
      tags$h4("1. Install CRAN packages"),
      
      tags$pre(
        "install.packages(c(
          'optparse',
          'jsonlite',
          'dplyr',
          'tibble',
          'stringr',
          'tidyr',
          'purrr',
          'Matrix',
          'data.table',
          'igraph',
          'FSelectorRcpp',
          'NetWeaver',
          'praznik',
          'reticulate',
          'openxlsx',
          'readxl',
          'Hmisc'
        ))"
      ),
      
      tags$h4("2. Install Bioconductor packages"),
      
      tags$pre(
        "if (!requireNamespace('BiocManager', quietly = TRUE))
          install.packages('BiocManager')"
      ),
      tags$pre(
        "BiocManager::install(c(
          'AnnotationDbi',
          'org.Hs.eg.db',
          'org.Mm.eg.db',
          'BiocParallel',
          'annotate'
        ))"
      ),
      
      tags$h4("Run command", style = "margin-top: 14px;"),
      tags$pre(
        style = "background:#0b1220; color:#e6edf3; padding:14px; border-radius:10px; font-size:13px;",
        "Rscript cli/DiCE_CLI.R --config cli/example_CLI.json"
      ),
      
      tags$hr(),
      
      # Inputs format
      tags$h3("Input file formats"),
      
      tags$h4("1) Differential gene expression analysis file (required)"),
      tags$ul(
        style = "font-size:15px; line-height:1.7;",
        tags$li("Supported: .csv / .tsv / .xlsx / .rds"),
        tags$li("Must include a gene column (e.g., Gene, Gene.Name, Gene.Symbol)"),
        tags$li("Must include: logFC, P.Value, adj.P.Val (column names should match your DiCE reader expectations)")
      ),
      
      tags$h4("2) Normalized gene expression matrix (required)"),
      tags$ul(
        style = "font-size:15px; line-height:1.7;",
        tags$li("Supported: .csv / .tsv / .xlsx / .rds"),
        tags$li("Rows = samples, Columns = genes"),
        tags$li("The last column must be the class/condition label (e.g., Tumor/Normal)")
      ),
      
      tags$div(
        class = "phase-summary",
        tags$p(tags$b("Tip:")),
        tags$p("Make sure treatment/control labels exactly match the values in the last column of the expression matrix (case-sensitive).")
      ),
      
      tags$hr(),
      
      # params
      h3("CLI Parameter Reference"),
      
      p(
        "The DiCE CLI uses a JSON configuration file to define analysis parameters. ",
        "The table below summarizes the supported parameters, their meaning, and default values."
      ),
      
      tags$table(
        class = "table table-bordered table-striped",
        style = "max-width:1100px;",
        
        tags$thead(
          tags$tr(
            tags$th("Parameter"),
            tags$th("Description"),
            tags$th("Allowed values / Example"),
            tags$th("Default")
          )
        ),
        
        tags$tbody(
          
          tags$tr(
            tags$td("species"),
            tags$td("Species used for gene annotation and STRING interaction mapping."),
            tags$td("'human' or 'mouse'"),
            tags$td("human")
          ),
          
          tags$tr(
            tags$td("dge_file_path"),
            tags$td("Path to the differential gene expression results file."),
            tags$td("Must contain Gene.Symbol, logFC, P.Value, adj.P.Val"),
            tags$td("Required")
          ),
          
          tags$tr(
            tags$td("normGeneExp_file_path"),
            tags$td("Path to the normalized gene expression matrix used for DiCE analysis."),
            tags$td("Genes × Samples. Genes under `Gene` column."),
            tags$td("Required")
          ),
          
          tags$tr(
            tags$td("metadata_file_path"),
            tags$td("Path to the metadata file containing sample annotations."),
            tags$td("Samples × Metadata columns. Include `sample_id` and `phenotype` columns."),
            tags$td("Required")
          ),
          
          tags$tr(
            tags$td("treatment"),
            tags$td("Label of treatment / case samples."),
            tags$td("Example: 'Tumor'"),
            tags$td("Required")
          ),
          
          tags$tr(
            tags$td("control"),
            tags$td("Label of control / reference samples."),
            tags$td("Example: 'Normal'"),
            tags$td("Required")
          ),
          
          tags$tr(
            tags$td("remove_pc_genes"),
            tags$td("Whether to remove non-protein-coding genes before downstream analysis."),
            tags$td("TRUE or FALSE."),
            tags$td("TRUE")
          ),
          
          tags$tr(
            tags$td("loose_criteria"),
            tags$td("Statistical significance metric used in Phase I filtering."),
            tags$td("'P.Value' or 'adj.P.Val'"),
            tags$td("adj.P.Val")
          ),
          
          tags$tr(
            tags$td("loose_cutoff"),
            tags$td("Threshold applied to the selected significance metric in Phase I."),
            tags$td("Numeric, e.g. 0.05"),
            tags$td("0.05")
          ),
          
          tags$tr(
            tags$td("logFC_cutoff"),
            tags$td("Minimum absolute log2 fold-change threshold for Phase I filtering."),
            tags$td("Numeric, e.g. 0"),
            tags$td("0")
          ),
          
          tags$tr(
            tags$td("ig_method"),
            tags$td("Method used for feature selection in Phase II."),
            tags$td("'IG', 'WIG', or 'none'"),
            tags$td("IG")
          ),
          
          tags$tr(
            tags$td("B"),
            tags$td("Number of bootstrap resamples used when ig_method = 'WIG'."),
            tags$td("Integer, e.g. 300"),
            tags$td("300")
          ),
          
          tags$tr(
            tags$td("ig_cutoff"),
            tags$td("Rule used to retain genes based on Information Gain."),
            tags$td("'all_mean', 'all_median', 'nonzero_mean', 'nonzero_median', 'all_nonzero', 'custom'"),
            tags$td("all_mean")
          ),
          
          tags$tr(
            tags$td("ig_custom_cutoff"),
            tags$td("Custom IG threshold used only when ig_cutoff = 'custom'."),
            tags$td("Numeric"),
            tags$td("NULL")
          ),
          
          tags$tr(
            tags$td("corr_method"),
            tags$td("Correlation method used in Phase III to calculate expression-based edge weights."),
            tags$td("'pearson' or 'spearman'"),
            tags$td("pearson")
          ),
          
          tags$tr(
            tags$td("ppi_db"),
            tags$td("Protein–protein interaction (PPI) database used to construct the network."),
            tags$td("Options: `stringdb` or `biogrid`."),
            tags$td("stringdb")
          ),
          
          tags$tr(
            tags$td("stringDB_confidence"),
            tags$td("Minimum STRINGdb confidence score threshold used to retain PPI interactions."),
            tags$td("Integer between 0–1000 (commonly 400 or 700). Used only when `ppi_db = 'stringdb'`."),
            tags$td("400 if ppi_db = `stringdb`")
          ),
          
          tags$tr(
            tags$td("centrality_list"),
            tags$td("List of centrality measures to compute in Phase IV."),
            tags$td("Example: ['betweenness', 'eigenvector']"),
            tags$td("betweenness, eigenvector")
          ),
          
          tags$tr(
            tags$td("dice_rules"),
            tags$td("A list of DiCE selection rules used to classify genes as final DiCE genes after ensemble ranking. Each rule can be created using helper functions such as dice_centrality_rule() or dice_ensemble_rule()."),
            tags$td("List"),
            tags$td(
              HTML("list(<br/>
               &nbsp;&nbsp;dice_centrality_rule(<br/>
               &nbsp;&nbsp;&nbsp;&nbsp;metric = \"betweenness\",<br/>
               &nbsp;&nbsp;&nbsp;&nbsp;threshold_type = \"percent\",<br/>
               &nbsp;&nbsp;&nbsp;&nbsp;threshold = 25<br/>
               &nbsp;&nbsp;),<br/>
               &nbsp;&nbsp;dice_centrality_rule(<br/>
               &nbsp;&nbsp;&nbsp;&nbsp;metric = \"eigen vector\",<br/>
               &nbsp;&nbsp;&nbsp;&nbsp;threshold_type = \"percent\",<br/>
               &nbsp;&nbsp;&nbsp;&nbsp;threshold = 25<br/>
               &nbsp;&nbsp;)<br/>
               )")
            )
          ),
          
          tags$tr(
            tags$td("dice_logic"),
            tags$td("Character string specifying how multiple dice_rules are combined to determine final DiCE genes."),
            tags$td("Character"),
            tags$td("AND")
          ),

          tags$tr(
            tags$td("modules"),
            tags$td("Configuration settings for module detection analysis on final DiCE genes."),
            tags$td("Nested JSON object containing module analysis parameters."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$run"),
            tags$td("Whether to run module detection analysis."),
            tags$td("TRUE or FALSE."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$unweighted"),
            tags$td("Settings for unweighted module detection."),
            tags$td("Nested JSON object."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$unweighted$enabled"),
            tags$td("Whether to run unweighted module detection."),
            tags$td("TRUE or FALSE."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$unweighted$algorithm"),
            tags$td("Community detection algorithm for unweighted modules."),
            tags$td("Options: `louvain` or `leiden`."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$unweighted$resolution"),
            tags$td("Resolution parameter controlling module granularity."),
            tags$td("Numeric value."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$unweighted$leiden_itrs"),
            tags$td("Number of Leiden optimization iterations."),
            tags$td("Integer value. Used only when algorithm = `leiden`."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$unweighted$leiden_beta"),
            tags$td("Beta parameter for Leiden optimization."),
            tags$td("Numeric value. Used only when algorithm = `leiden`."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$weighted"),
            tags$td("Settings for weighted module detection."),
            tags$td("Nested JSON object."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$weighted$enabled"),
            tags$td("Whether to run weighted module detection."),
            tags$td("TRUE or FALSE."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$weighted$algorithm"),
            tags$td("Community detection algorithm for weighted modules."),
            tags$td("Options: `louvain` or `leiden`."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$weighted$resolution"),
            tags$td("Resolution parameter controlling module granularity."),
            tags$td("Numeric value."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$weighted$leiden_itrs"),
            tags$td("Number of Leiden optimization iterations."),
            tags$td("Integer value. Used only when algorithm = `leiden`."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("modules$weighted$leiden_beta"),
            tags$td("Beta parameter for Leiden optimization."),
            tags$td("Numeric value. Used only when algorithm = `leiden`."),
            tags$td("Optional")
          ),
          
          tags$tr(
            tags$td("seed"),
            tags$td("Random seed for reproducibility."),
            tags$td("Integer, e.g. 123"),
            tags$td("123")
          ),
          
          tags$tr(
            tags$td("outdir"),
            tags$td("Output directory where CLI results and logs are saved."),
            tags$td("Example: 'dice_cli_output'"),
            tags$td("dice_cli_output")
          ),
          
          tags$tr(
            tags$td("job_name"),
            tags$td("Prefix used for naming output files."),
            tags$td("Example: 'DiCE_CLI_Run'"),
            tags$td("DiCE_CLI_Run")
          )
          
        )
      ),
      br(),
      
      # Config JSON
      tags$h3("Example config JSON"),
      
      tags$p(
        style = "font-size:15px; color:#555;",
        "Save as ",
        tags$code("cli/example_CLI.json"),
        " (update paths to your local machine)."
      ),
      
      tags$pre(
        style = "background:#0b1220; color:#e6edf3; padding:14px; border-radius:10px; font-size:13px; white-space:pre-wrap;",
        paste0(
          '{
            "species": "human",
            "dge_file_path": "sample_data/sample_Human_data_DGE.Rds",
            "normGeneExp_file_path": "sample_data/sample_Human_data_geneExp.Rds",
            "metadata_file_path": "sample_data/sample_Human_metadata.csv",
            "treatment": "Tumor",
            "control": "Normal",
            "loose_criteria": "adj.P.Val",
            "loose_cutoff": 0.05,
            "logFC_cutoff": 0,
            "ig_method": "IG",
            "B": 300,
            "ig_cutoff": "all_mean",
            "ig_custom_cutoff": null,
            "corr_method": "pearson",
            "ppi_db": "stringdb",
            "stringDB_confidence": 400,
            "remove_pc_genes": false,
            "centrality_list": "betweenness,eigenvector",
            "dice_rules": [
              {
                "type": "centrality",
                "metric": "betweenness",
                "threshold_type": "percent",
                "threshold": 25
              },
              {
                "type": "ensemble",
                "threshold_type": "rank",
                "threshold": 200
              }
            ],
            "dice_logic": "OR",
            "modules": {
              "run": true,
              
              "unweighted": {
                "enabled": true,
                "algorithm": "louvain",
                "resolution": 1,
                "leiden_itrs": 10,
                "leiden_beta": 0.01
              },
              
              "weighted": {
                "enabled": true,
                "algorithm": "louvain",
                "resolution": 1,
                "leiden_itrs": 10,
                "leiden_beta": 0.01
              }
            },
            "seed": 123,
            "outdir": "dice_cli_output",
            "job_name": "DiCE_CLI_Run"
          }'

        )
      ),
      
      tags$hr(),
      
      # Outputs
      tags$h3("Outputs"),
      
      tags$ul(
        style = "font-size:15px; line-height:1.7;",
        tags$li(tags$b("DiCE results Excel:"), " ", tags$code("<outdir>/<job_name>_dice_results.xlsx")),
        tags$li(tags$b("Modules Excel (if enabled):"), " ", tags$code("<outdir>/<job_name>_modules.xlsx"))
      ),
      
      tags$hr(),
      
      tags$p(
        tags$a(
          href = "#",
          onclick = "Shiny.setInputValue('go_about_tab', Math.random())",
          style = "font-weight:600; color:#2c7fb8;",
          "See more details on the Home page →"
        )
      ),
      tags$hr(),
      
      # Troubleshooting
      tags$h3("Troubleshooting"),
      
      tags$h4("Common issues"),
      tags$ul(
        style = "font-size:15px; line-height:1.7;",
        tags$li(tags$b("Parameters are NULL:"), " usually JSON key names don't match the CLI option names."),
        tags$li(tags$b("File not found:"), " paths in JSON must be absolute or correct relative paths."),
        tags$li(tags$b("Group labels mismatch:"), " treatment/control must match class labels in expression matrix exactly."),
        tags$li(tags$b("Missing packages:"), " install required packages in the same R environment used by Rscript.")
      ),
      
      tags$h4("Recommended run (prints to terminal + writes log)"),
      tags$pre(
        style = "background:#0b1220; color:#e6edf3; padding:14px; border-radius:10px; font-size:13px;",
        "Rscript cli/DiCE_CLI.R --config cli/example_CLI.json 2>&1 | tee dice_cli_output/DiCE_CLI_Run_terminal.log"
      ),
      
      tags$div(style = "height: 24px;")
    )
  )
}