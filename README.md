# Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry (CyTOF) data

This repository contains R scripts to reproduce the analyses and figures in our paper comparing clustering methods for high-dimensional flow cytometry and mass cytometry (CyTOF) data.

- Weber L.M. and Robinson M.D. (2016) *Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data.* Cytometry Part A, 89A: 1084–1096. [Open access available here](http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.23030/full).

    Note: Supplementary files (Supporting Information PDF and Supporting Information Table S1) are available via the journal link above, or from the [supplementary_files](supplementary_files/) directory (latest version: November 18, 2016).



## Overview

In this study, we compared the performance of 18 clustering methods for automated detection of cell populations in high-dimensional flow cytometry and mass cytometry (CyTOF) data, using 6 publicly available data sets from experiments in immunology as benchmarks. These results extend previously published comparisons by focusing on high-dimensional data and including new methods developed for CyTOF data.

A list of the clustering methods can be found in Table 1 in the paper. A list of the data sets is provided in Table 2.

The comparisons showed that several methods performed well, including FlowSOM, X-shift, PhenoGraph, Rclusterpp, and flowMeans. Among these, FlowSOM had extremely fast runtimes, making this method well-suited for interactive, exploratory analysis of large, high-dimensional data sets on a standard laptop or desktop computer.

Based on our results, we recommend the use of FlowSOM (with manual selection of the number of clusters; see paper) as a first choice for this type of analysis, since this method gave best or near-best performance across all data sets, together with extremely fast runtimes.

See the paper for more details, in particular regarding the advantages of the different methods for different clustering tasks (detecting multiple cell populations vs. detecting a single rare population).



## FlowSOM and Rtsne example code

FlowSOM [(Van Gassen et al., 2015)](http://www.ncbi.nlm.nih.gov/pubmed/25573116) is available as a [Bioconductor package](http://bioconductor.org/packages/release/bioc/html/FlowSOM.html) for the R statistical programming language.

A worked example showing how to use FlowSOM for clustering and [Rtsne](https://github.com/jkrijthe/Rtsne) for visualization is available in the [FlowSOM-Rtsne-example](https://github.com/lmweber/FlowSOM-Rtsne-example) repository.



## Updates

Updated results for new clustering algorithms or new reference data sets will be published on this website.

The following updates are currently available:

- densityCut: New clustering method [(Ding et al., 2016)](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btw227). R scripts and a summary report are included in the [updates/densityCut](updates/densityCut/) directory.



## Contents of this repository

R scripts are organized into the following subdirectories. All scripts contain comments explaining the overall purpose and individual steps.

- [data_preparation](data_preparation/): preparation of benchmark data files
- [ensemble_clustering](ensemble_clustering/): run and evaluate ensemble clustering
- [evaluate_results](evaluate_results/): scripts to evaluate results from all methods
- [helpers](helpers/): helper functions
- [plots_and_tables](plots_and_tables/): generate plots and tables of results
- [range_k](range_k/): run and evaluate FlowSOM over range of values k (number of clusters)
- [run_methods](run_methods/): scripts to run all methods (or instructions to run graphical interfaces, where required)
- [stability_analysis](stability_analysis/): run and evaluate methods for stability analysis

Supplementary files from the published paper are included in the following directory:

- [supplementary_files](supplementary_files/): supplementary files from paper (latest version: November 18, 2016)

R scripts and summary reports for updated results are included in the following directory:

- [updates](updates/): updated results for new clustering methods or new reference data sets



## Data files

Pre-processed data files for the benchmark data sets are available from FlowRepository ([repository FR-FCM-ZZPH](https://flowrepository.org/id/FR-FCM-ZZPH)).

Original data files can be obtained through the references listed in Table 2 in the paper.



