# Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry (CyTOF) data

This repository contains R scripts to reproduce the analyses and figures in our paper on comparing clustering methods for high-dimensional flow cytometry and mass cytometry (CyTOF) data.

A preprint of the paper is available on [bioRxiv](http://biorxiv.org/content/early/2016/04/07/047613).


## Overview

In this study, we compared the performance of 18 clustering methods for automated detection of cell populations in high-dimensional flow cytometry and mass cytometry (CyTOF) data sets. A list of the clustering methods can be found in Table 1 in the paper.

The comparisons showed that FlowSOM (with manual selection of the number of clusters) gave the best performance across data sets, and had extremely fast runtimes. We recommend FlowSOM as a first choice for this type of analysis. In particular, the fast runtimes enable interactive, exploratory analyses of large data sets (possibly up to millions of cells) on a standard laptop or desktop computer.

Several other methods also performed well across multiple data sets, including X-shift, PhenoGraph, Rclusterpp, and flowMeans.



## FlowSOM and Rtsne example code

FlowSOM [(van Gassen et al. 2016)](http://www.ncbi.nlm.nih.gov/pubmed/25573116) is available as a [Bioconductor package](http://bioconductor.org/packages/release/bioc/html/FlowSOM.html) for the R programming language. A worked example showing how to use FlowSOM for clustering and [Rtsne](https://github.com/jkrijthe/Rtsne) for visualization is available in the [FlowSOM-Rtsne-example](https://github.com/lmweber/FlowSOM-Rtsne-example) repository.



## Contents

Scripts are organized into the following subfolders. All scripts contain comments explaining their purpose and the individual steps.

- [data_preparation](data_preparation/): preparation of benchmark data files
- [ensemble_clustering](ensemble_clustering/): run and evaluate ensemble clustering
- [evaluate_results](evaluate_results/): scripts to evaluate results from all methods
- [helpers](helpers/): helper functions
- [plots_and_tables](plots_and_tables/): generate plots and tables of resutls
- [range_k](range_k/): run and evaluate methods over range of values for k (number of clusters)
- [run_methods](run_methods/): scripts to run all methods (or instructions to run graphical interfaces, where required)
- [stability_analysis](stability_analysis/): run and evaluate methods for stability analysis



## Data files

Pre-processed data files for each of the benchmark data sets are available from FlowRepository ([repository FR-FCM-ZZPH](https://flowrepository.org/id/FR-FCM-ZZPH)).

Original data files can be obtained through the references in Table 2 in the paper.


