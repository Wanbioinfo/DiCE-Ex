# DiCE
## Version: V1.2.3

## Disease Biomarker Identification with scRNA-seq Data

### Required R Libraries

```r
dplyr, stringr, tibble, tidyr, purrr, Matrix, parallel, data.table, igraph, FSelectorRcpp, NetWeaver, NewWave, zinbwave, SingleCellExperiment, 
SummarizedExperiment, AnnotationDbi, annotate, org.Hs.eg.db, org.Mm.eg.db, BiocParallel, stats, utils, openxlsx, readxl, praznik, scuttle, Hmisc
```
---


### Installation

```r
install.packages("path/to/DiCE_1.2.3.tar.gz", repos = NULL, type = "source")

```
---

### Step 1: Prepare Input Files for DiCE

If you are running DiCE on **bulk RNA-seq** data, you need the following two input files:

1. **"dge_data"** – Rds/xlsx/csv/tsv file containing the differential gene expression analysis results of **protein coding genes**.
  - The columns must include: **"Gene.Symbol", "logFC", "P.Value", "adj.P.Val"**. 
  - Note: Column names in your data must exactly match the above.

2. **"logNorm_geneExp"** - Rds/xlsx/csv/tsv file containing the normalized (logCPM) gene expression data of **protein coding genes**.
  - Rows represent samples, and columns represent genes.
  - There should be one additional column indicating the sample group/class.
    - Example: If a sample is from a tumor tissue, the Group column should contain “Tumor”.


If you are running DiCE on **scRNA-seq** data, you need the following three input files:


1. **"dge_data"** – Rds/xlsx/csv/tsv file containing the differential gene expression analysis results of **protein coding genes**.
  - The columns must include: **"Gene.Symbol", "logFC", "P.Value", "adj.P.Val"**. 
  - Note: Column names in your data must exactly match the above.

2. **"logNorm_geneExp"** – Rds/xlsx/csv/tsv file containing the normalized (logCPM) gene expression data of **protein coding genes**.
  - Rows represent samples, and columns represent genes.
  - There should be one additional column indicating the sample group/class.
    - Example: If a sample is from a tumor tissue, the Group column should contain “Tumor”.

3. **"raw_gene_exp"** – Rds/xlsx/csv/tsv file containing the raw UMI counts of the **protein coding genes** .
  - Rows represent cells, and columns represent genes.
  - There should be one additional column indicating the cell group/class.
    - Example: If a sample is from a tumor tissue, the Group column should contain “Tumor”.

---

### Step 2: Run `DiCE.R` to Identify DiCE Genes

Configure the parameters based on your application:

- **data_type** : Sequencing data type ("bulkRNA-seq" OR "scRNA-seq"). Default is "bulkRNA-seq".
- **species** : "human" or "mouse". Default is "human".
- **dge_file_path** : File path to the differential gene expression file. The columns must include: "Gene.Symbol", "logFC", "P.Value", "adj.P.Val".
- **normGeneExp_file_path** : File path for the differential gene expression analysis '.Rds' file. (cells/samples x genes + label column). When the data_type = "scRNA-seq" and ig_method = "none", you do not need to input this file.
- **rawGeneExp_file_path** : File path for the raw UMI counts file (Needs only when the data_type = "scRNA-seq"). (cells/samples x genes + label column)
- **treatment** : Label of treatment samples (eg: Tumor)
- **control** : Label of control samples (eg: Normal)
- **remove_pc_genes** : If TRUE, protein-coding genes are removed from the input data before analysis. Default TRUE.
- **loose_criteria** : Statistical significance metric used for initial gene filtering ("P.Value" or "adj.P.Val"). Default is "adj.P.Val"
- **loose_cutoff** : Numeric threshold for filtering based on 'loose_criteria'. Genes with values ≤ this cutoff are retained. Default is 0.05.
- **logFC_cutoff** : Minimum absolute log2 fold change threshold for retaining genes. Default is 0.
- **ig_method** : Character string specifying the Information Gain strategy to use. Options are:
	- "IG" : Standard information gain (Default)
	- "wIG" : Weighted information gain
	- "none" : Skip information gain filtering
- **B** : Number of bootstrap resamples used in weighted IG calculation. Default is 300.
- **ig_cutoff** - Method for selecting IG-filtered genes. Options include: 
	- "all_mean" : Retain all genes with IG greater than the mean IG computed across all genes (including zeros). Default.
	- "all_median" : Retain all genes with IG greater than the median IG computed across all genes (including zeros)
	- "nonzero_mean" : Compute the mean IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero mean threshold
	- "nonzero_median" : Compute the median IG only among genes with IG > 0, and retain genes whose IG exceeds this non-zero median threshold
	- "all_nonzero" : Retain all genes with IG > 0
	- "custom" : Retain genes with IG greater than a user-defined threshold specified in **ig_custom_cutoff**
- **ig_custom_cutoff** : Numeric IG threshold used only when **ig_cutoff = "custom"**.
- **corr_mode** : Mode for computing gene–gene correlation. Determines how gene expression data are pre-processed and how correlation is computed. Options:
	- "directCorr" : Use raw or normalized expression values without dropping zero-expression cells. Computes correlation directly. (Default)
	- "remove_Zerocells" : Exclude cell pairs where both genes have zero expression before computing correlation. Useful for sparse data. (For scRNAseq data).
	- "ZINB-WaVE" : Apply ZINB-WaVE denoising to model zero inflation and overdispersion before computing correlation. (For scRNAseq data). 
- **corr_method** : "pearson" OR "spearman". Default is "pearson".
- **corr_pval_cutoff** : P-value threshold (default 0.05); correlations with p-value > cutoff are set to 0.
- **centrality_list** : Character vector of centrality metrics to compute. Valid options include: "betweenness", "eigen vector", "pagerank", "closeness", "harmonic", "authority", "strength". 
-**dice_rules** : List of rules used to identify final DiCE genes after centrality calculation and ensemble ranking. Rules should be created using helper functions such as 					dice_centrality_rule() and dice_ensemble_rule(). For a centrality rule, a gene passes if it satisfies the specified cutoff in either the treatment or control network 						for that centrality. For an ensemble rule, a gene passes if its Ensemble_Rank satisfies the specified cutoff. Default uses Betweenness top 25% in either treatment or 						control network.
-**dice_logic** : Logic used to combine multiple rules in dice_rules. Options include:
					"AND" : A gene must pass all rules to be classified as a DiCE gene. Default.
					"OR" : A gene must pass at least one rule to be classified as a DiCE gene.
- **markTF**: If TRUE, transcription factor genes will be marked. Default FALSE.

---
### Instructions to define DiCE rules

**dice_rules** are used to define how final DiCE genes are selected after centrality calculation and ensemble ranking. Rules can be created using `dice_centrality_rule()` for centrality-based filtering and `dice_ensemble_rule()` for ensemble-rank-based filtering.

Each rule is evaluated separately. 

A gene is considered to pass a centrality-based rule if it satisfies the cutoff in either the treatment network or the control network for the selected centrality metric. For example, if the rule is based on Betweenness, the logic is:

`(Betweenness_treatment passes cutoff) | (Betweenness_control passes cutoff)`

For `dice_centrality_rule(metric, threshold_type, threshold)`:
- **metric** specifies the centrality measure to evaluate, such as "betweenness", "eigenvector", "pagerank", or other supported centrality metrics.
- **threshold_type** specifies how the cutoff should be interpreted. The cutoff depends on threshold_type:
		- "mean": the gene must have a centrality value greater than or equal to the mean centrality value in either treatment or control. In this case, **threshold** is not required.
		- "percent": the gene must have a centrality value greater than or equal to the cutoff defined by the top specified percentage in either treatment or control. For example, 					threshold = 25 means the gene must be greater than or equal to the 75th percentile in either network.
		- "rank": the gene must fall within the top specified rank in either treatment or control. Here, **threshold** should be a rank value such as 200.
- **threshold** provides the numeric cutoff value when required.

An `dice_ensemble_rule(threshold_type, threshold)` rule is evaluated using Ensemble_Rank. A gene passes if its ensemble **rank** is within the specified cutoff.

When multiple rules are provided, they are combined using **dice_logic**:
- "AND": the gene must pass all rules
- "OR": the gene must pass at least one rule

A gene that satisfies the final combined logic is labeled as DiCE; otherwise, it keeps its previous phase label.

#### Example 1: single centrality rule

```r
dice_rules <- list(
  dice_centrality_rule(
    metric = "betweenness",
    threshold_type = "mean"
  )
)
dice_logic <- "AND"
```
This means a gene is classified as DiCE if:
`(Betweenness_treatment >= mean(Betweenness_treatment)) | (Betweenness_control  >= mean(Betweenness_control))`

#### Example 2: multiple centrality rules
```r
dice_rules <- list(
  dice_centrality_rule("betweenness", "mean"),
  dice_centrality_rule("eigenvector", "mean")
)
dice_logic <- "AND"
```
This means a gene is classified as DiCE if:

`((Betweenness_treatment >= mean(Betweenness_treatment)) | (Betweenness_control  >= mean(Betweenness_control))) & ((EigenVector_treatment >= mean(EigenVector_treatment)) | (EigenVector_control  >= mean(EigenVector_control)))`

#### Example 3: centrality rule plus ensemble rule

```r
dice_rules <- list(
  dice_centrality_rule("betweenness", "percent", 25),
  dice_ensemble_rule("rank", 200)
)
dice_logic <- "AND"
```
This means a gene is classified as DiCE if it is:
	- within the top 25% of Betweenness in either treatment or control, and
	- within the top 200 genes by Ensemble_Rank

#### Default rule

If no custom rule is provided, the default DiCE rule is:

```r
dice_rules <- list(
  dice_centrality_rule(
    metric = "betweenness",
    threshold_type = "percent",
    threshold = 25
  ),
  dice_centrality_rule(
    metric = "eigen vector",
    threshold_type = "percent",
    threshold = 25
  )
)
dice_logic <- "AND"
```

---
### DiCE Output

- **dice_results_df**: A data frame containing genes evaluated across Phase I, Phase II, Phase III, and final DiCE selection, including differential expression statistics, centrality measures, and final ensemble rankings.

- **interactions_df**: A data frame of protein–protein interaction pairs used for network construction, including their corresponding edge weights (absolute correlation coefficients).

- **phase2_denoised_geneExp (scRNA-seq only)**: A filtered denoised gene expression matrix used for Phase II network construction and correlation analysis.

---

### DiCE usage examples

Exmaple: Run DiCE pipeline on bulk RNA-seq data

```r
dice_rules = list(
  dice_centrality_rule(
    metric = "betweenness",
    threshold_type = "percent",
    threshold = 25
    
  ),
  dice_centrality_rule(
    metric = "eigen vector",
    threshold_type = "percent",
    threshold = 25
  )
)
dice_logic <- "AND"

dice_results <- perform_DiCE(
  data_type = "bulkRNA-seq",
  species = "human",
  dge_file_path = "path/to/dge_results.rds",
  normGeneExp_file_path = "path/to/logNorm_geneExp.rds",
  treatment = "Tumor",
  control = "Normal",
  loose_criteria = "adj.P.Val",
  loose_cutoff = 0.05,
  logFC_cutoff = 1,
  ig_method = "IG",
  corr_mode = "directCorr",
  corr_method = "pearson",
  dice_rules = dice_rules,
  dice_logic = dice_logic,
  markTF = TRUE
)

# View DiCE genes
head(dice_results$dice_results_df)
```

Exmaple: Run DiCE pipeline on scRNA-seq data

```r
dice_rules = list(
  dice_centrality_rule(
    metric = "betweenness",
    threshold_type = "percent",
    threshold = 25
    
  ),
  dice_centrality_rule(
    metric = "eigen vector",
    threshold_type = "percent",
    threshold = 25
  )
)
dice_logic <- "AND"

dice_results <- perform_DiCE(
  data_type = "scRNA-seq",
  species = "human",
  dge_file_path = "path/to/dge_results.rds",
  normGeneExp_file_path = "path/to/logNorm_geneExp.rds",
  rawGeneExp_file_path = "path/to/raw_geneExp.rds",
  treatment = "Tumor",
  control = "Normal",
  loose_criteria = "adj.P.Val",
  loose_cutoff = 0.05,
  logFC_cutoff = 1,
  ig_method = "IG",
  corr_mode = "ZINB-WaVE",
  corr_method = "spearman",
  dice_rules = dice_rules,
  dice_logic = dice_logic,
  markTF = TRUE
)

# View DiCE genes
head(dice_results$dice_results_df)
```

**Note**: For scRNA-seq data, if are not using Information Gain (ig_method = none), you do not need to input normalized gene expression data. Providing differential gene expression data and gene raw counts is enough. 

---

### (Optional): Identify **condition-sepcific weighted** PPI network modules,

`detect_PPI_weightedModules()` identifies condition-specific modules in weighted PPI networks for a given gene set. It constructs treatment and control networks using condition-specific edge weights, detects communities, and summarizes module structure and node-level roles.

- **Input:** Gene list, PPI interaction table with condition-specific edge weights, and names of treatment and control conditions.
	- For the PPI interaction table, you can use `dice_results$interactions_df`. 
	- Otherwise, the table must contain the columns `source`, `target`, `weight_Tumor`, and `weight_Normal` when the treatment is `Tumor` and the control is `Normal`.
- **Output:** Module summaries (`network_modules`), module statistics (`module_stats`), module overlap between conditions (Jaccard similarity `cmp_jaccard` and common gene counts `cmp_count`), node-level module role statistics (`nodes_stats`), all the edges with modules information(`edges_with_modules`), and condition specific intermodule connectivity
(`treatment_interMod`, `control_interMod`).
	

**Node-level module role statistics**  (`nodes_stats`) table describes the network role of each gene within the detected modules.

	- **Condition** – Condition in which the network role was calculated (e.g., Tumor or Normal). 
	- **Gene.Symbol** – Gene corresponding to the network node. 
	- **Module** – Module (community) the gene belongs to in that condition. 
	- **Within_Module_Z** – Within-module degree z-score. Measures how strongly the gene is connected to other genes within its own module. Higher values indicate hub-like genes within 								the module.
	- **Participation_Coefficient** – Measures how evenly a gene’s connections are distributed across different modules. Higher values indicate connector genes linking multiple modules.
	- **Intra_Module_Degree** – Weighted degree of connections within the gene’s own module.
	- **Inter_Module_Degree** – Weighted degree of connections from the gene to genes in other modules.


Configure the parameters for `detect_PPI_weightedModules()`:

- **gene_list** : Character vector of genes to retain in the networks.
- **interactions_df** : Data frame of interactions with `source`, `target`, and condition-specific weight columns named `weight_<condition>` (e.g., `weight_Tumor`).
- **treatment** : Character string naming the treatment condition (used to select `weight_<treatment>`. E.g., "Tumor").
- **control** : Character string naming the control condition (used to select `weight_<control>`. E.g., "Normal").
- **louvain_resolution** : Resolution parameter for computing modularity in Louvain algorithm (default 1).
- **seed** : Random seed for reproducibility (default 123).


Example to find condition-sepcific weighted modules.

```r
res <- detect_PPI_weightedModules(
			gene_list = my_genes,
			interactions_df = intr_df,
			treatment = "Tumor",
			control = "Normal",
			louvain_resolution = 0.9,
			seed = 123
)

# View network stats of the nodes in each module
head(res$nodes_stats)

```

### (Optional): Identify **unweighted** PPI network modules (not a condition-specific analysis)

`detect_PPI_unweightedModules()` identifies community modules in an unweighted STRING v12 PPI network for a given gene set. It extracts STRING interactions for the selected genes (combined score ≥ 400), builds the PPI graph, detects Louvain communities, and returns module summaries, gene-to-module assignments, and intra-module interaction tables.

- **Input:**
  - `gene_list`: genes for the module analysis
  - `species`: `"human"` or `"mouse"` (loads the corresponding STRING v12 reference files)
  - `seed`: random seed for reproducibility

- **Output:** 
  - `summary_df`: number of modules and modularity
  - `membership_df`: gene-to-module assignments and within-module degree
  - `edges_by_module`: intra-module interaction tables split by module

Example to find unweighted network modules.

```r
# Run module detection on DiCE genes
modules <- detect_PPI_unweightedModules(
     gene_list = my_genes,
     species = "human",
     seed = 123
)

# View module membership
head(modules$membership_df)

```
### Updates

#### Last update - 16/03/2026

#### Version 1.2.3
- 16/03/2026 : Added flexible rule-based final DiCE gene selection using dice_rules and dice_logic, with support for custom centrality/ensemble cutoffs and rule-tracking columns in the results table.

#### Version 1.2.3
- 09/03/2026 : Included P-value threshold (default 0.05); correlations with p-value > cutoff are set to 0 (Used Hmisc rcorr()).

#### Version 1.2.1
- 06/03/2026 - Condition-specific modules returns a dataframe with edges with correspoding modules and weights, and summary of condition-specific intermodule connectivity.
- 06/03/2026 - Temporarily removed NewWave,

#### Version 1.2.0
- 04/03/2026 -  Added a function to find condition-specific modules. 
- 27/02/2026 -  Added option to annotate the transcription factor genes
- 27/02/2026 -  Change the script so that normalized gene expression data is not required for scRNA-seq data when ig_method = "none".
- 25/02/2026 -  Replaced ZINB-WaVE with NewWave for single-cell RNA-seq preprocessing, improving computational efficiency and scalability.

#### Version 1.1.5
- 25/02/2026 -  Added changes such that DiCE can skip Phase2
- 20/02/2026 -  Set IG cutoff value option to filter genes in Phase2

#### Version 1.1.4
- 11/02/2026 - Removing non-protein coding genes moved as an option to choose.
- 05/12/2025 - DiCE now supports input files in multiple formats, including .rds, .xlsx, .tsv, and .csv, allowing users to load data seamlessly regardless of file type. 
- 05/12/2025 - DiCE now supports input files in multiple formats, including .rds, .xlsx, .tsv, and .csv, allowing users to load data seamlessly regardless of file type. 
- 13/11/2025 - Added different IG cutoffs (all_mean, all_median, nonzero_mean, nonzero_median, all_nonzero)
- 13/11/2025 - Added a function at the beginning of the DiCE pipeline to retain only protein-coding genes before downstream analysis.
- 03/11/2025 - Introduced other network centralities (Authority, Strength, Closeness, pagerank, Harmonic). User can select the centralities they want.

#### Version 1.1.3
- 21/10/2025 - Included function to find the network modules/communities among the DiCE genes with unweighted PPI and Louvain algorithm. 
- 15/10/2025 - Introduced class weighted Information Gain with Monte carlo approximation.

#### Version 1.1.2
- 12/10/2025 - Handled errors/warnings due to column names mismatches in the input data.
- 10/10/2025 - Forced correlations which were undefined due to zero variance to be zero.

#### Version 1.1.1
- 08/10/2025 - Added IG information for Phase 1 genes, updated the script to return only the dataframe containing all Phase 1, Phase 2, Phase 3, and DiCE genes with their final ensemble ranking, and handled variations in column names provided by the user.

#### Version 1.1.0
- 20/09/2025 - Corrected the final ensemble ranking calculation of all genes.
- 17/09/2025 - Changed the PPI construction method to use downloaded StringDB files instead of the StringDB R package and creates the final gene ranking table. 
- 10/09/2025 - Corrected eigen vector centrality calculation using correlation coefficient as the edge weights.

			
