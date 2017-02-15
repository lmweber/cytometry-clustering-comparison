#########################################################################################
# R script to load and evaluate results from all methods (including densityCut)
#
# Lukas Weber, February 2017
#########################################################################################


# note: some of these scripts take a few minutes to run


DIR_PREVIOUS <- "../../../evaluate_results"
DIR_DENSITYCUT <- getwd()


source("evaluate_densityCut.R")


setwd(DIR_PREVIOUS)


source("evaluate_ACCENSE.R")
source("evaluate_ClusterX.R")
source("evaluate_DensVM.R")
source("evaluate_FLOCK.R")
source("evaluate_flowClust.R")
source("evaluate_flowMeans.R")
source("evaluate_flowMerge.R")
source("evaluate_flowPeaks.R")
source("evaluate_FlowSOM.R")
source("evaluate_FlowSOM_pre.R")
source("evaluate_immunoClust.R")
source("evaluate_kmeans.R")
source("evaluate_PhenoGraph.R")
source("evaluate_Rclusterpp.R")
source("evaluate_SamSPECTRAL.R")
source("evaluate_SPADE.R")
source("evaluate_SWIFT.R")
source("evaluate_Xshift.R")


setwd(DIR_DENSITYCUT)


# update data directory (for densityCut)

DATA_DIR <- "../../../../../benchmark_data_sets"



