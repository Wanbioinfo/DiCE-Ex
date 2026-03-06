about_tab <- function() {
  tabPanel(
    title = "About",
    value = "about",
    div(
      class = "page-container",
      
      h2("About DiCE"),
      p(
        "DiCE (Differential Centrality-Ensemble analysis) is a multi-phase, network-guided gene discovery framework designed to identify 
        influential disease-associated genes that may not exhibit large fold-changes in traditional differential expression analysis. 
        By integrating Information Gain (IG), condition-specific weighted PPI networks, and multiple network centrality measures, 
        DiCE highlights genes whose functional influence shifts between biological states. This approach has been validated in 
        cancer datasets and successfully identifies survival-linked and cancer-fitness genes overlooked by classical differential gene expression (DGE) methods."
      ),
      
      br(),
      
      div(
        class = "phase-summary",
        
        tags$h4("Publication"),
        
        tags$p(
          tags$b("DiCE: differential centrality-ensemble analysis based on gene expression profiles and protein-protein interaction network")
        ),
        
        tags$p(
          "Pashaei, E., Liu, S., Li, K.Y., Zang, Y., Yang, L., Lautenschlaeger, T., Huang, J., Lu, X., & Wan, J. (2025). DiCE: differential centrality-ensemble analysis based on gene expression profiles and protein–protein interaction network. Nucleic Acids Research, 53."
        ),
        
        tags$p(
          tags$a(
            "View article (Nucleic Acids Research)",
            href   = "https://academic.oup.com/nar/article/53/13/gkaf609/8192812",
            target = "_blank"
          )
        )
      ),
      
      br(),
      
      tags$div(
        class = "phase-summary",
        tags$h4("Bug reports and feature requests"),
        tags$p(
          style = "font-size:14px;",
          "Please use our GitHub issue tracker to report bugs, request new features, or ask questions about DiCE-Ex."
        ),
        tags$a(
          href   = "https://github.com/Wanbioinfo/DiCE_Ex/issues",
          target = "_blank",
          class  = "btn btn-primary",
          style  = "font-size:14px; margin-top:10px;",
          "Report an Issue on GitHub"
        )
      ),
      
      br(),
      
      div(
        id = "docs_section",
        
        h2("Documentation"),
        
        h3("Overview"),
        p(
          "DiCE-Ex implements the complete DiCE workflow, allowing users to upload bulk RNA-seq datasets and 
          run all phases—from DGE-based candidate gene selection to IG filtering, network reconstruction, centrality analysis, 
          and ensemble ranking. The platform supports human and mouse datasets and provides interactive visualization of 
          DiCE genes, PPI subnetworks, and community modules."
        ),
        p(
          "DiCE-Ex provides clear advantages over standalone R packages by offering an integrated, end-to-end, and user-configurable workflow."
        ),
        
        br(),
        
        h3("Input Requirements"),
        p("DiCE-Ex supports multiple file formats including .csv, .tsv, .xlsx, and .rds."),
        
        h4("Bulk RNA-seq Inputs"),
        
        tags$ul(
          tags$li(
            tags$p(
              "Normalized expression matrix (logCPM), formatted as samples × genes with one additional class/Group/Condition column. ",
              "Do not place sample names as the first column."
            ),
            tags$img(
              src = "images/gene_expr_image.png",
              style = "max-width:600px; margin-top:6px; border:1px solid #ddd; border-radius:6px;"
            )
          ),
          tags$li(
            tags$p("Differential expression analysis (DEA) results with columns: Gene.Symbol, P.Value, adj.P.Val, logFC."),
            tags$img(
              src = "images/dge_image.png",
              style = "max-width:300px; margin-top:6px; border:1px solid #ddd; border-radius:6px;"
            )
          )
        ),
        
        br(),
        
        h3("Workflow Summary"),
        p("DiCE consists of five main phases that progressively refine and prioritize biologically meaningful genes:"),
        
        h4("Phase I — Significance Filtering"),
        
        p(
          "Phase I constructs the initial candidate gene pool using differential gene expression (DGE) results supplied by the user. 
          This step applies loose statistical filters to retain genes that may show modest but biologically meaningful expression changes."
        ),
        
        tags$ul(
          
          tags$li(
            tags$b("P-value measure: "),
            "Users can select which statistical measure should be used for filtering. 
            The default option is ",
            tags$b("Adjusted p-value"),
            ", which accounts for multiple testing correction. 
            Alternatively, the raw ",
            tags$b("P.Value"),
            " column can be used if adjusted values are not available."
          ),
          
          tags$li(
            tags$b("Threshold for selected p-value: "),
            "Defines the significance cutoff applied to the selected p-value measure. 
            Genes with values below this threshold are retained in the candidate pool. 
            The default value (0.05) is intentionally permissive to capture genes with subtle expression differences."
          ),
          
          tags$li(
            tags$b("Threshold for |log2FC|: "),
            "Defines the minimum absolute log2 fold-change required for a gene to pass Phase I filtering. 
            A value of ",
            tags$b("0"),
            " means that genes can pass the filter even if they show very small fold changes, 
            allowing DiCE to detect genes whose importance arises from network topology rather than large expression shifts."
          )
        ),
        
        h4("Phase II — Information Gain (IG) Filtering"),
        
        p(
          "Phase II identifies genes whose expression profiles best distinguish the two biological conditions. 
          This step applies Information Gain (IG), an entropy-based feature selection metric commonly used in 
          machine learning. IG measures how much information a gene provides about class labels, allowing DiCE 
          to retain genes with strong discriminative power even if their fold-change is modest."
        ),
        
        tags$ul(
          
          tags$li(
            tags$b("Information Gain method: "),
            "Users can select between ",
            tags$b("Information Gain (IG)"),
            ", ",
            tags$b("Weighted Information Gain (WIG)"),
            ", or ",
            tags$b("None"),
            ". The standard IG method evaluates how well each gene separates the two conditions. 
            The WIG option applies Monte Carlo resampling to balance class sizes before computing IG, 
            which improves robustness when sample groups are imbalanced. Selecting ",
            tags$b("None"),
            " skips Phase II filtering and retains all Phase I genes."
          ),
          
          tags$li(
            tags$b("B (bootstrap resamples for WIG): "),
            "Defines the number of bootstrap resampling iterations used when computing Weighted Information Gain. 
            During each iteration, samples from the larger class are randomly subsampled to match the size of 
            the smaller class, and IG is recomputed. The final WIG score is obtained by averaging across all 
            bootstrap iterations. The default value (300) provides stable estimates while maintaining 
            reasonable computational efficiency."
          ),
          
          tags$li(
            tags$b("IG cutoff rule: "),
            "Defines how genes are filtered based on their IG scores. Several automatic rules are available:",
            
            tags$ul(
              tags$li(tags$b("IG > mean (all genes): "), 
                      "Retains genes whose IG score is greater than the average IG across all genes."),
              
              tags$li(tags$b("IG > median (all genes): "), 
                      "Uses the median IG value across all genes as the cutoff."),
              
              tags$li(tags$b("IG > mean (IG > 0): "), 
                      "Computes the mean IG value using only genes with non-zero IG scores."),
              
              tags$li(tags$b("IG > median (IG > 0): "), 
                      "Computes the median IG value using only genes with non-zero IG scores."),
              
              tags$li(tags$b("IG > 0: "), 
                      "Retains all genes with non-zero IG values."),
              
              tags$li(tags$b("Custom: "), 
                      "Allows users to define their own IG threshold.")
            )
          )
        ),
        
        p(
          tags$b("Note: "),
          "Information Gain is calculated using the gene expression matrix and the phenotype labels. 
          Genes that pass the IG filtering step proceed to the network-based analysis in Phase III."
        ),
        
        h4("Phase III — PPI Network Construction"),
        
        p(
          "In Phase III, DiCE reconstructs condition-specific protein–protein interaction (PPI) networks 
          using the candidate genes retained from Phase II. These networks integrate known biological 
          interactions from the STRING database V12.0 with expression-derived correlation weights, enabling 
          detection of phenotype-specific network rewiring."
        ),
        
        tags$ul(
          
          tags$li(
            tags$b("STRING interaction network: "),
            "DiCE uses curated protein–protein interactions from the STRING database (Version 12.0). 
            Only interactions between genes present in the candidate gene set are retained. 
            Genes that do not appear in STRING or have no interactions within the candidate set 
            are removed from the network."
          ),
          
          tags$li(
            tags$b("Condition-specific networks: "),
            "Separate PPI networks are constructed for each biological condition (e.g., treatment 
            and control) using the corresponding subset of samples from the expression matrix."
          ),
          
          tags$li(
            tags$b("Edge weights: "),
            "Interaction edges are weighted using expression correlations between genes. ",
            "For each gene pair, the correlation coefficient is computed using samples belonging to a specific condition. ",
            "The edge weight is defined as ",
            tags$b("|correlation|"),
            ", while the distance between two genes is defined as ",
            tags$b("1 - |correlation|"),
            "."
          ),
          
          tags$li(
            tags$b("Correlation type: "),
            "Users can choose between ",
            tags$b("Pearson"),
            " and ",
            tags$b("Spearman"),
            " correlation when calculating edge weights. ",
            
            tags$ul(
              tags$li(
                tags$b("Pearson correlation: "),
                "Measures linear relationships between gene expression profiles and is 
                appropriate when expression values follow approximately normal distributions."
              ),
              
              tags$li(
                tags$b("Spearman correlation: "),
                "Measures rank-based associations and is more robust to outliers and 
                non-linear relationships."
              )
            )
          )
          
        ),
        
        p(
          tags$b("Note: "),
          "By constructing separate weighted PPI networks for each condition, DiCE enables 
          comparison of network topology between biological states, which forms the basis 
          for the centrality-based prioritization performed in Phase IV."
        ),
        
        h4("Phase IV — Network Centrality Analysis"),
        tags$ul(
          tags$li(
            "Computes multiple centrality measures selected from betweenness, eigenvector, closeness, PageRank, harmonic, strength, and authority.",
            
            tags$div(
              style = "margin-top:10px; margin-left:15px;",
              
              tags$p(
                tags$b("Note: "),
                "Avoid selecting highly correlated centrality measures together, as they capture similar topological information."
              ),
              
              tags$p(
                "For example, Eigenvector centrality and Authority often exhibit strong correlations, 
                indicating redundant topological information."
              ),
              
              fluidRow(
                column(
                  width = 6,
                  tags$img(
                    src   = "images/tcga_prad_heatmap.png",
                    style = "width:100%; border:1px solid #ddd; border-radius:6px;",
                    alt   = "Centrality correlations in Cancer Genome Atlas Prostate Adenocarcinoma (TCGA-PRAD)"
                  ),
                  tags$p(
                    style = "text-align:center; font-size:0.9em; margin-top:5px;",
                    "TCGA-PRAD dataset"
                  )
                ),
                column(
                  width = 6,
                  tags$img(
                    src   = "images/nepc_heatmap.png",
                    style = "width:100%; border:1px solid #ddd; border-radius:6px;",
                    alt   = "Centrality correlations in a Neuroendocrine Prostate Cancer (NEPC) dataset"
                  ),
                  tags$p(
                    style = "text-align:center; font-size:0.9em; margin-top:5px;",
                    "NEPC 2019 dataset"
                  )
                )
              )
            )
          ),
          
          tags$li(
            "Genes that exceed user-defined thresholds (mean, median, or top-K%) in at least one condition for all selected centralities are designated as ",
            tags$b("DiCE genes.")
          )
        ),
        
        h4("Phase V — Ensemble Ranking"),
        tags$ul(
          tags$li("All genes are ranked based on the absolute differences in each selected centrality measure between the two conditions."),
          tags$li("A consensus ensemble rank is then computed using the ProductOfRank method.")
        ),
        
        br(),
        
        h3("Outputs and Result Interpretation"),
        
        tags$p(
          "DiCE-Ex outputs are designed to support both gene-level prioritization and network-level interpretation. 
          Below we describe how users should interpret each type of result returned by the platform."
        ),
        
        h4("DiCE Results Table Columns"),
        
        p(
          "The DiCE results table integrates differential expression statistics, feature selection scores, 
          network centrality measures, and ensemble rankings. The following table describes each column."
        ),
        
        tags$table(
          class = "table table-bordered table-striped",
          style = "max-width:1000px;",
          
          tags$thead(
            tags$tr(
              tags$th("Column"),
              tags$th("Description")
            )
          ),
          
          tags$tbody(
            
            tags$tr(
              tags$td("Gene.Symbol"),
              tags$td("Official gene symbol used as the primary identifier in the analysis.")
            ),
            
            tags$tr(
              tags$td("logFC"),
              tags$td("Log2 fold change between treatment and control groups obtained from differential expression analysis.")
            ),
            
            tags$tr(
              tags$td("P.Value"),
              tags$td("Raw p-value from the differential expression analysis.")
            ),
            
            tags$tr(
              tags$td("adj.P.Val"),
              tags$td("Adjusted p-value corrected for multiple testing (e.g., Benjamini–Hochberg FDR).")
            ),
            
            tags$tr(
              tags$td("IG / Weighted IG"),
              tags$td("Information Gain score measuring how well a gene separates the two biological conditions.")
            ),
            
            tags$tr(
              tags$td("Centrality measures for each condition"),
              tags$td("Network topology metrics computed for each gene in the condition-specific PPI networks. Examples include betweenness, eigenvector, PageRank, closeness, harmonic, strength, and authority.")
            ),
            
            tags$tr(
              tags$td("Centrality differences between the two conditions"),
              tags$td("Absolute differences in centrality values between treatment and control networks, capturing changes in gene influence between conditions.")
            ),
            
            tags$tr(
              tags$td("Centrality_rank columns"),
              tags$td("Rank of the gene based on the absolute difference in a specific centrality measure between treatment and control networks. ",
                      "Genes with larger centrality differences receive higher priority (lower rank values).")
            ),
            
            tags$tr(
              tags$td("pass_count"),
              tags$td("Number of selected centrality measures for which the gene exceeds the user-defined filtering threshold ",
                          "(e.g., mean, median, or top-K%).")
            ),
            
            tags$tr(
              tags$td("pass_centralities"),
              tags$td(
                "Lists the centrality measures for which the gene satisfied the DiCE filtering criteria. ",
                "For example, 'Betweenness, EigenVector' indicates the gene passed the threshold for both metrics."
              )
            ),
            
            tags$tr(
              tags$td("Phase"),
              tags$td("Indicates the stage of the DiCE pipeline in which the gene was retained. Genes labeled as 'DiCE' represent the final prioritized gene set.")
            ),
            
            tags$tr(
              tags$td("ProductOfRank"),
              tags$td(
                "Intermediate score computed by multiplying the ranks of the gene across all selected centrality measures. ",
                "Genes consistently ranked highly across multiple centralities obtain smaller ProductOfRank values."
              )
            ),
            tags$tr(
              tags$td("Ensemble_Rank"),
              tags$td(
                "Final consensus ranking of genes derived from the ProductOfRank score. ",
                "Genes with lower Ensemble_Rank values represent stronger candidate genes with consistent network influence changes."
              )
            )
            
          ),
          tags$img(
            src = "images/dice_results.png",
            style = "max-width:1000px; width:100%; margin-top:15px; border:1px solid #ddd; border-radius:6px;"
          ),
          
          tags$p(
            style = "font-size:0.9em; color:#555; text-align:center; margin-top:5px;",
            "Example DiCE results table showing gene ranking and filtering options. 
            Users can search for genes, filter by phase, and inspect network-based prioritization metrics."
          )
        ),
        
        h5("How to Interpret DiCE Gene Rankings"),
        
        tags$ul(
          
          tags$li(
            "Genes ranked at the top of the table typically show both expression relevance and strong changes in network influence between conditions."
          ),
          
          tags$li(
            "Some prioritized genes may not exhibit large fold-changes but demonstrate significant shifts in network centrality, highlighting genes that may regulate condition-specific network rewiring."
          ),
          
          tags$li(
            tags$b("Phase column interpretation: "),
            "The Phase column indicates the stage of the DiCE workflow in which the gene satisfied filtering criteria. ",
            "Genes labeled as ",
            tags$b("DiCE"),
            " passed the centrality-based filtering step (Phase IV) and represent the final DiCE gene set."
          ),
          
          tags$li(
            "Genes labeled as Phase III were retained during earlier filtering steps but did not satisfy all DiCE centrality criteria. ",
            "These genes may still appear highly ranked because the final ranking considers centrality differences across all genes."
          ),
          
          tags$li(
            "Therefore, some Phase III genes may appear above DiCE genes in the ranking table. ",
            "This does not indicate a contradiction: the Phase column reflects filtering criteria, whereas the ranking reflects the magnitude of network topology changes."
          ),
          
          tags$li(
            "Users are encouraged to investigate both top-ranked DiCE genes and highly ranked Phase III genes, particularly when exploring pathway enrichment or network modules."
          )
          
        ),
        
        h3("Interpretation of DiCE Modules"),
        
        p(
          "DiCE-Ex identifies community modules within the STRING-derived protein–protein interaction (PPI) network 
          constructed from DiCE genes. Modules represent groups of genes that interact more frequently with 
          each other than with the rest of the network and may correspond to shared biological pathways or processes."
        ),
        
        tags$ul(
          
          tags$li(
            tags$b("Module summary: "),
            "Displays the number of detected modules and the modularity score, which measures how strongly the network is divided into communities."
          ),
          
          tags$img(
            src = "images/module_summary.png",
            style = "max-width:1000px; width:50%; margin-top:15px; border:1px solid #ddd; border-radius:6px;"
          ),
          
          tags$p(
            style="font-size:0.9em; color:#555; text-align:left;",
            "Example module summary table. In this application, 5 modules were detected and the modularity is 0.3258."
          ),
          
          tags$li(
            tags$b("Module membership table: "),
            "Lists genes belonging to each module along with network connectivity metrics such as within-module degree."
          ),
          
          tags$img(
            src = "images/module_memb.png",
            style = "max-width:1000px; width:50%; margin-top:15px; border:1px solid #ddd; border-radius:6px;"
          ),
          
          tags$p(
            style="font-size:0.9em; color:#555; text-align:left;",
            "Example gene membership table in Module 2. The table shows the genes and their degree in the Module 2."
          ),
          
          tags$li(
            tags$b("Network visualization: "),
            "Shows module-specific subnetworks where nodes represent genes and edges represent protein–protein interactions. 
            Selecting a module or clicking on a gene displays the corresponding visualization."
          ),
          
          tags$img(
            src = "images/net_vis.png",
            style = "max-width:1000px; width:50%; margin-top:15px; border:1px solid #ddd; border-radius:6px;"
          ),
          
          tags$p(
            style="font-size:0.9em; color:#555; text-align:left;",
            "DiCE PPI module network of gene FAM72B. Nodes represent genes and edges represent protein–protein interactions."
          ),
          
          tags$li(
            tags$b("Gene exploration: "),
            "Clicking a gene in the network displays expression boxplots comparing treatment and control conditions."
          ),
          tags$img(
            src = "images/gene_explore.png",
            style = "max-width:1000px; width:50%; margin-top:15px; border:1px solid #ddd; border-radius:6px;"
          ),
          
          tags$p(
            style="font-size:0.9em; color:#555; text-align:left;",
            "Gene exploration panel for the gene FAM72B."
          ),
          
        ),
        
        br(),
        
        tags$ul(
          tags$li("Genes within the same module often participate in related biological functions or pathways."),
          tags$li("Users can perform Gene Ontology (GO) or pathway enrichment analysis on genes within each module ",
                  "to identify enriched biological processes and better understand their functional relevance.")
        ),
        
        br(),
        
      ),
      
      br(),
      div(
        id = "support_section",
        h2("Support"),
        tags$ul(
          tags$li(
            tags$b("Report bugs / request features: "),
            tags$a(
              href = "https://github.com/Wanbioinfo/DiCE_Ex/issues",
              target = "_blank",
              "GitHub Issues"
            )
          ),
          tags$li(
            tags$b("Please include: "),
            "job title, dataset name (if public), a screenshot of the error, and the log output from the Results page."
          )
        )
      ),
      br(),
      
    ) # end page-container
  )
}
