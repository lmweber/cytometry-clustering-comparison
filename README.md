# Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data

**NOTE: We are currently working on an updated version of this paper. R scripts in this repository have only been partially updated so far (as at Aug 29th, 2016). The updated version should be complete by around early to mid September, so please check back soon.**

This repository contains R scripts to reproduce the analyses and figures in our paper on comparing clustering methods for high-dimensional single-cell flow cytometry and mass cytometry (CyTOF) data.

A preprint of the paper is available on [bioRxiv](http://biorxiv.org/content/early/2016/04/07/047613).


## Overview

In this study, we compared the performance of 13 clustering methods (11 distinct methods and two variations) for automatically detecting clusters representing cell populations in high-dimensional flow cytometry and mass cytometry (CyTOF) data sets. A list of the clustering methods can be found in Table 1 in the paper.

We evaluated the clustering methods using four publicly available data sets from experiments in immunology, with cell population identities available from manual gating. A summary of the data sets is provided in Table 2 in the paper.

The comparisons showed that FlowSOM (available as a [Bioconductor package](http://bioconductor.org/packages/release/bioc/html/FlowSOM.html); link to [paper](http://www.ncbi.nlm.nih.gov/pubmed/25573116)) performed well across all data sets, and had extremely fast runtimes. We recommend FlowSOM as a first choice for analyzing new data sets from experiments in high-dimensional flow and mass cytometry.


## Contents

All R scripts are saved in the [scripts](scripts/) folder.

Each script contains comments explaining the purpose of the script and the individual steps within it.

The main R script to reproduce the figures shown in the paper is [scripts/plots_main_results.R](scripts/plots_main_results.R). R scripts to reproduce the data preprocessing steps are saved in [scripts/data_preparation_scripts](scripts/data_preparation_scripts/).


## Data files

Preprocessed data files required to reproduce the analyses are available from FlowRepository ([repository FR-FCM-ZZPH](https://flowrepository.org/id/FR-FCM-ZZPH)).

Original data files are available from the references listed in Table 2 in the paper.


