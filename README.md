# Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data

This repository contains R scripts to reproduce the analysis and figures in our paper on comparing clustering methods for high-dimensional single-cell flow and mass cytometry data.

(reference and link to bioRxiv preprint here when submitted)


## Overview

In this study, we compared the performance of 13 clustering methods (11 distinct methods and two variations) for automatically detecting clusters representing cell populations in high-dimensional flow and mass cytometry data sets. A list of the clustering methods can be found in Table 1 in the paper.

We tested the clustering methods using four publicly available, real data sets from high-dimensional flow and mass cytometry experiments in immunology. A summary of the data sets is provided in Table 2 in the paper.

(sentence on main results)


## Contents

All R scripts are saved in the [scripts](scripts/) folder.

Each script contains comments explaining the overall purpose of the script, as well as the individual steps within it.

The main R script to reproduce the figures shown in the paper is [scripts/plots_main_results.R](scripts/plots_main_results.R). R scripts to reproduce the data preprocessing steps are saved in [scripts/data_preparation_scripts](scripts/data_preparation_scripts/).


## Data files

Preprocessed data files required to reproduce the analysis and figures are available from FlowRepository (repository FR-FCM-ZZPH).

Original data files are publicly available from the references listed in Table 2 in the paper.


