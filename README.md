# Project Overview

The iDMET project provides a framework for interpreting metabolomic data using differential metabolite profiles. This repository extends iDMET to perform enrichment analysis on metabolite sets.

The main goals are:

1. Process and organize metabolomic datasets from multiple studies.
2. Identify significant metabolites in each comparison.
3. Perform odds ratio calculations and enrichment analysis.
4. Generate annotation lists for MSEA (Metabolite Set Enrichment Analysis).
5. Visualize metabolite networks and pathway enrichment results.

## Requirement

- igraph
- cytoscape
- clusterProfiler

## Dataset Description
* Annoation list (metabodic.csv)

  * The annotation list presents metabolites in an organized manner, listing each common metabolite name together with its synonyms. This structured information facilitates enrichment analysis and helps ensure consistent identification of metabolites across different datasets.

* Data
  
  * The data used for iDMET analysis, consisting of 130 studies, can be found in the `data` folder.  
    The data used as queries for the enrichment analysis are located in the `query_data` folder.

## References
- Rira Matsuta, Hiroyuki Yamamoto, Masaru Tomita, Rintaro Saito, "iDMET: network-based approach for integrating differential analysis of cancer metabolomics", BMC Bioinformatics, 23:508 (2022).


