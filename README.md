# DiCE-Ex
**DiCE-Ex** is an interactive web server implementing the **Differential Centrality-Ensemble (DiCE)** framework for network-based gene discovery from transcriptomic data. The platform enables researchers to identify influential disease-associated genes by integrating differential expression, information-theoretic feature selection, and protein–protein interaction (PPI) network topology—without relying on prior disease annotations.

- Web server link: [https://wanbioinfo-dice-ex.share.connect.posit.cloud/](https://connect.posit.iu.edu/DiCE-Ex/)
---

## Background

Traditional differential gene expression (DGE) analysis often misses biologically important genes with subtle expression changes. Network-based approaches help address this limitation but frequently depend on prior disease knowledge.

**DiCE (Differential Centrality-Ensemble analysis)** overcomes these challenges by:
- Integrating expression-driven gene selection with network topology
- Constructing condition-specific weighted PPI networks
- Ranking genes based on changes in network influence

DiCE-Ex provides a complete, user-friendly interface to perform this analysis directly from transcriptomic data.

---

## Key Features

- **Multi-phase DiCE pipeline**
  - Differential expression filtering
  - Information Gain (IG) or Weighted IG–based feature selection
  - STRING-based PPI network construction
  - Differential network centrality analysis

- **Interactive analysis**
  - Ranked gene tables with filtering and search
  - PPI module detection
  - Interactive PPI subnetworks and gene-level expression visualization

- **Asynchronous job execution**
  - Real-time progress tracking
  - Shareable, bookmarkable result links
  - No login required

---

## Input Requirements

### 1. Normalized Expression Matrix
- File formats: `.csv`, `.tsv`, `.xlsx`, `.rds`
- Genes as **rows**
- Samples as **columns**
- **Do not keep genes as rownames. Keep them under `Gene` column.**
- Do **not** place sample names as the first column

### 2. Differential Gene Expression (DGE) Results
Must include the following columns:
- Gene identifier (official gene symbols recommended)
- `logFC`
- `P.Value`
- `adj.P.Val`

### 3. Metadata File
- File formats: `.csv`, `.tsv`, `.xlsx`, `.rds`
- Samples as **rows**
- Must contain:
  - `sample_id` column matching sample names in the expression matrix
  - `phenotype` column containing group/class labels
- Group labels are case-sensitive

### 4. Experimental Design
- Two-group comparison (e.g., Tumor vs Normal)
- Treatment and control labels must exactly match values in the `phenotype` column of the metadata file (case-sensitive)

---

## Supported Species

- **Human (Homo sapiens)**
- **Mouse (Mus musculus)**

---

## Output

- DiCE results containing all genes retained across DiCE phases, along with differential expression metrics, IG or Weighted IG scores, network centrality values for each condition, and final ensemble ranks.
- Detected PPI modules and membership tables: Detected modules in the STRING-derived interaction network of DiCE genes represent groups of tightly connected genes that may correspond to shared biological processes or pathways.
- Interactive results visualization and exploration
- Downloadable result files (Excel)

---

## Publication

Pashaei, E., Liu, S., Li, K.Y., Zang, Y., Yang, L., Lautenschlaeger, T., Huang, J., Lu, X., & Wan, J. (2025). DiCE: differential centrality-ensemble analysis based on gene expression profiles and protein–protein interaction network. Nucleic Acids Research, 53. 

https://academic.oup.com/nar/article/53/13/gkaf609/8192812

## DiCE R package

https://github.com/Wanbioinfo/DiCE_package

