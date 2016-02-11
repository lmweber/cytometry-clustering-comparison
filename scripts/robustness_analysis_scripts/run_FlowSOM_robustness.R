#########################################################################################
# R script to run FlowSOM for robustness analysis
# 
# This script runs one iteration of FlowSOM for the robustness analysis. The main script 
# "robustness_analysis.R" then runs the scripts for each clustering method several times
# in a loop.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the robustness analysis.
# 
# Lukas M. Weber, February 2016
#########################################################################################


library(flowCore)
library(FlowSOM)

# helper functions
source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for robustness analysis
source("load_results_truth_robustness.R")



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")


# FlowSOM requires input data as flowFrame objects

data_Levine_32 <- flowCore::read.FCS(file_Levine_32, transformation = FALSE)
data_Levine_13 <- flowCore::read.FCS(file_Levine_13, transformation = FALSE)
data_Nilsson <- flowCore::read.FCS(file_Nilsson, transformation = FALSE)
data_Mosmann <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)

# head(data_Levine_32)
# head(data_Levine_13)
# head(data_Nilsson)
# head(data_Mosmann)

# dim(data_Levine_32)
# dim(data_Levine_13)
# dim(data_Nilsson)
# dim(data_Mosmann)


# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- 7:21

# length(marker_cols_Levine_32)
# length(marker_cols_Levine_13)
# length(marker_cols_Nilsson)
# length(marker_cols_Mosmann)




###################
### Run FlowSOM ###
###################

# run FlowSOM

# set.seed(123)
runtime_FlowSOM_Levine_32 <- system.time({
  fSOM_Levine_32 <- FlowSOM::ReadInput(data_Levine_32, transform = FALSE, scale = FALSE)
  fSOM_Levine_32 <- FlowSOM::BuildSOM(fSOM_Levine_32, colsToUse = marker_cols_Levine_32)
  fSOM_Levine_32 <- FlowSOM::BuildMST(fSOM_Levine_32)
})


# set.seed(123)
runtime_FlowSOM_Levine_13 <- system.time({
  fSOM_Levine_13 <- FlowSOM::ReadInput(data_Levine_13, transform = FALSE, scale = FALSE)
  fSOM_Levine_13 <- FlowSOM::BuildSOM(fSOM_Levine_13, colsToUse = marker_cols_Levine_13)
  fSOM_Levine_13 <- FlowSOM::BuildMST(fSOM_Levine_13)
})


# set.seed(123)
runtime_FlowSOM_Nilsson <- system.time({
  fSOM_Nilsson <- FlowSOM::ReadInput(data_Nilsson, transform = FALSE, scale = FALSE)
  fSOM_Nilsson <- FlowSOM::BuildSOM(fSOM_Nilsson, colsToUse = marker_cols_Nilsson)
  fSOM_Nilsson <- FlowSOM::BuildMST(fSOM_Nilsson)
})


# set.seed(123)
runtime_FlowSOM_Mosmann <- system.time({
  fSOM_Mosmann <- FlowSOM::ReadInput(data_Mosmann, transform = FALSE, scale = FALSE)
  fSOM_Mosmann <- FlowSOM::BuildSOM(fSOM_Mosmann, colsToUse = marker_cols_Mosmann)
  fSOM_Mosmann <- FlowSOM::BuildMST(fSOM_Mosmann)
})


# plots

# FlowSOM::PlotStars(fSOM_Levine_32)
# FlowSOM::PlotStars(fSOM_Levine_13)
# FlowSOM::PlotStars(fSOM_Nilsson)
# FlowSOM::PlotStars(fSOM_Mosmann)


# extract cluster labels

# str(fSOM_Levine_32$map)

# head(fSOM_Levine_32$map$mapping)
# dim(fSOM_Levine_32$map$mapping)


clus_FlowSOM_Levine_32 <- fSOM_Levine_32$map$mapping[, 1]
clus_FlowSOM_Levine_13 <- fSOM_Levine_13$map$mapping[, 1]
clus_FlowSOM_Nilsson <- fSOM_Nilsson$map$mapping[, 1]
clus_FlowSOM_Mosmann <- fSOM_Mosmann$map$mapping[, 1]

# length(clus_FlowSOM_Levine_32)
# length(clus_FlowSOM_Levine_13)
# length(clus_FlowSOM_Nilsson)
# length(clus_FlowSOM_Mosmann)


# cluster sizes and number of clusters

# table(clus_FlowSOM_Levine_32)
# table(clus_FlowSOM_Levine_13)
# table(clus_FlowSOM_Nilsson)
# table(clus_FlowSOM_Mosmann)

# length(table(clus_FlowSOM_Levine_32))
# length(table(clus_FlowSOM_Levine_13))
# length(table(clus_FlowSOM_Nilsson))
# length(table(clus_FlowSOM_Mosmann))




########################
### Run FlowSOM_meta ###
########################

# run optional metaclustering step (FlowSOM_meta)

# set number of clusters
# k_Levine_32 <- 20
# k_Levine_13 <- 30
# k_Nilsson <- 50
# k_Mosmann <- 50


# set.seed(123)
# runtime_FlowSOM_meta_Levine_32 <- system.time({
#   meta_clustering_Levine_32 <- FlowSOM::metaClustering_consensus(fSOM_Levine_32$map$codes, k = k_Levine_32)
# })


# set.seed(123)
# runtime_FlowSOM_meta_Levine_13 <- system.time({
#   meta_clustering_Levine_13 <- FlowSOM::metaClustering_consensus(fSOM_Levine_13$map$codes, k = k_Levine_13)
# })


# set.seed(123)
# runtime_FlowSOM_meta_Nilsson <- system.time({
#   meta_clustering_Nilsson <- FlowSOM::metaClustering_consensus(fSOM_Nilsson$map$codes, k = k_Nilsson)
# })


# set.seed(123)
# runtime_FlowSOM_meta_Mosmann <- system.time({
#   meta_clustering_Mosmann <- FlowSOM::metaClustering_consensus(fSOM_Mosmann$map$codes, k = k_Mosmann)
# })


# alternatively: set number of clusters automatically (does not perform well)

# set.seed(123)
# runtime_FlowSOM_meta_Levine_32 <- system.time({
# meta_clustering_Levine_32 <- FlowSOM::MetaClustering(fSOM_Levine_32$map$codes, 
#                                                      method = "metaClustering_consensus", 
#                                                      max = 50)
# })
# 
# set.seed(123)
# runtime_FlowSOM_meta_Levine_13 <- system.time({
# meta_clustering_Levine_13 <- FlowSOM::MetaClustering(fSOM_Levine_13$map$codes, 
#                                                      method = "metaClustering_consensus", 
#                                                      max = 50)
# })
# 
# set.seed(123)
# runtime_FlowSOM_meta_Nilsson <- system.time({
#   meta_clustering_Nilsson <- FlowSOM::MetaClustering(fSOM_Nilsson$map$codes, 
#                                                     method = "metaClustering_consensus", 
#                                                     max = 50)
# })
# 
# set.seed(123)
# runtime_FlowSOM_meta_Mosmann <- system.time({
#   meta_clustering_Mosmann <- FlowSOM::MetaClustering(fSOM_Mosmann$map$codes, 
#                                                     method = "metaClustering_consensus", 
#                                                     max = 50)
# })


# combine runtime

# runtime_FlowSOM_meta_Levine_32 <- runtime_FlowSOM_Levine_32 + runtime_FlowSOM_meta_Levine_32
# runtime_FlowSOM_meta_Levine_13 <- runtime_FlowSOM_Levine_13 + runtime_FlowSOM_meta_Levine_13
# runtime_FlowSOM_meta_Nilsson <- runtime_FlowSOM_Nilsson + runtime_FlowSOM_meta_Nilsson
# runtime_FlowSOM_meta_Mosmann <- runtime_FlowSOM_Mosmann + runtime_FlowSOM_meta_Mosmann


# extract cluster labels

# meta_clustering_Levine_32
# meta_clustering_Levine_13
# meta_clustering_Nilsson
# meta_clustering_Mosmann

# clus_FlowSOM_meta_Levine_32 <- meta_clustering_Levine_32[fSOM_Levine_32$map$mapping[, 1]]
# clus_FlowSOM_meta_Levine_13 <- meta_clustering_Levine_13[fSOM_Levine_13$map$mapping[, 1]]
# clus_FlowSOM_meta_Nilsson <- meta_clustering_Nilsson[fSOM_Nilsson$map$mapping[, 1]]
# clus_FlowSOM_meta_Mosmann <- meta_clustering_Mosmann[fSOM_Mosmann$map$mapping[, 1]]

# length(clus_FlowSOM_meta_Levine_32)
# length(clus_FlowSOM_meta_Levine_13)
# length(clus_FlowSOM_meta_Nilsson)
# length(clus_FlowSOM_meta_Mosmann)


# cluster sizes and number of clusters

# table(clus_FlowSOM_meta_Levine_32)
# table(clus_FlowSOM_meta_Levine_13)
# table(clus_FlowSOM_meta_Nilsson)
# table(clus_FlowSOM_meta_Mosmann)

# length(table(clus_FlowSOM_meta_Levine_32))
# length(table(clus_FlowSOM_meta_Levine_13))
# length(table(clus_FlowSOM_meta_Nilsson))
# length(table(clus_FlowSOM_meta_Mosmann))




#########################
### CALCULATE RESULTS ###
#########################

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FlowSOM_Levine_32 <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine_32, clus_truth_Levine_32)
res_FlowSOM_Levine_13 <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine_13, clus_truth_Levine_13)
res_FlowSOM_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Nilsson, clus_truth_Nilsson)
res_FlowSOM_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Mosmann, clus_truth_Mosmann)

# res_FlowSOM_Levine_32
# res_FlowSOM_Levine_13
# res_FlowSOM_Nilsson
# res_FlowSOM_Mosmann




###############################################################
### OUTPUT RESULTS FOR ONE ITERATION OF ROBUSTNESS ANALYSIS ###
###############################################################

# data sets with multiple populations (Levine_32, Levine_13)
# output mean F1 score, mean precision, mean recall, runtime

res_robust_FlowSOM_Levine_32 <- list(
  mean_F1 = mean(res_FlowSOM_Levine_32$F1), 
  mean_pr = mean(res_FlowSOM_Levine_32$pr), 
  mean_re = mean(res_FlowSOM_Levine_32$re), 
  runtime = unname(runtime_FlowSOM_Levine_32["elapsed"])
)

res_robust_FlowSOM_Levine_13 <- list(
  mean_F1 = mean(res_FlowSOM_Levine_13$F1), 
  mean_pr = mean(res_FlowSOM_Levine_13$pr), 
  mean_re = mean(res_FlowSOM_Levine_13$re), 
  runtime = unname(runtime_FlowSOM_Levine_13["elapsed"])
)


# data sets with a single rare population of interest (Nilsson, Mosmann)
# output F1 score, precision, recall (for population of interest), and runtime

res_robust_FlowSOM_Nilsson <- list(
  F1 = as.numeric(res_FlowSOM_Nilsson$F1), 
  pr = as.numeric(res_FlowSOM_Nilsson$pr), 
  re = as.numeric(res_FlowSOM_Nilsson$re), 
  runtime = unname(runtime_FlowSOM_Nilsson["elapsed"])
)

res_robust_FlowSOM_Mosmann <- list(
  F1 = as.numeric(res_FlowSOM_Mosmann$F1), 
  pr = as.numeric(res_FlowSOM_Mosmann$pr), 
  re = as.numeric(res_FlowSOM_Mosmann$re), 
  runtime = unname(runtime_FlowSOM_Mosmann["elapsed"])
)


