# Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data

## Overview

This repository contains R scripts to reproduce the analysis in our paper on comparing clustering methods for high-dimensional single-cell flow and mass cytometry data.

We compared 13 clustering methods (including two variations of methods), using four publicly available, real data sets from experiments in immunology. The data sets are:

- Levine_2015_marrow_32: a 32-dimensional mass cytometry data set containing 14 major immune cell populations, from Levine et al. (2015).

- Levine_2015_marrow_13: a 13-dimensional mass cytometry data set containing 24 major immune cell populations, from Levine et al. (2015).

- Nilsson_2013_HSC: a 13-dimensional flow cytometry data set, containing a single rare cell population of interest, from Nilsson et al. (2013). The rare population consists of hematopoietic stem cells (HSCs), and accounts for 0.8% of total cells in the data set.

- Mosmann_2014_activ: a 15-dimensional flow cytometry data set, containing a single rare cell population of interest, from Mosmann et al. (2014). The rare population consists of activated (cytokine-producing) memory CD4 T cells, and accounts for 0.03% of total cells in the data set.


## Contents

R scripts are saved in the [scripts](scripts/) folder.

R scripts to prepare the benchmark data sets are saved in [scripts/data_preparation_scripts](scripts/data_preparation_scripts).

Currently, the repository contains R scripts only. The final version will also include PDF files of all figures included in the paper, as well as additional information in the README.
