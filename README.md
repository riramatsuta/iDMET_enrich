# iDMET+ Enrichment Analysis Pipeline

## Overview
The iDMET project provides a framework for interpreting metabolomic data using differential metabolite profiles.  
Steps 1–4 in the analysis workflow are preparatory steps used to generate the metabolite sets for enrichment analysis. 
To perform enrichment analysis with your own data, you only need to execute Step 5 (iDMET_msea_source.R, msea_ora.R).


The analysis consists of the following steps:

1. Process and organize metabolomic datasets from multiple studies  (standardizes metabolite names and prepares datasets):  `db_step1.R`
2. Identify significant metabolites in each comparison:  `db_step2.R`
3. Perform odds ratio calculations:  `db_step3.R`
4. Generate metabolite set lists:  `iDMET_setlist.R`
5. Enrichment analysis:  `iDMET_msea_source.R, msea_ora.R`

## Requirements
- R 4.0 or higher
- R packages:
  - igraph
  - cytoscape
  - clusterProfiler

## Dataset Description
* Annoation list (metabodic.csv)

  * The annotation list presents metabolites in an organized manner, listing each common metabolite name together with its synonyms. This structured information facilitates enrichment analysis and helps ensure consistent identification of metabolites across different datasets.

* Pathway Set List

  * The file `iDMET_pathwaylist.csv` contains the main pathway set list provided by this study.

- **Format**: CSV file
- **Columns**:
  - Pathway ID
  - Pathway Name
  - Metabolites included in the pathway
- **Usage**: This file is used by `msea_ora.R` to perform pathway-based enrichment analysis.  

* Data
  
  * The data used for iDMET analysis, consisting of 130 studies, can be found in the `data` folder.  
    The data used as queries for the enrichment analysis are located in the `query_data` folder.


##  How to Perform Enrichment Analysis
Use the script `msea_ora.R` to perform the enrichment analysis.  
Specify the paths to your CSV file and the output directory in the script.  

* Input Data
  * You need to prepare a CSV file containing **all metabolites detected in your dataset**, including both those that changed (differential metabolites) and those that did not. The file should include:
     * Metabolite names
     * Fold changes or p-values


This ensures that the enrichment analysis accounts for all detected metabolites, reducing potential bias that could arise from using only significantly changed metabolites.


* Output
  * Enrichment results (CSV)


## References
- Rira Matsuta, Hiroyuki Yamamoto, Masaru Tomita, Rintaro Saito, "iDMET: network-based approach for integrating differential analysis of cancer metabolomics", BMC Bioinformatics, 23:508 (2022).




