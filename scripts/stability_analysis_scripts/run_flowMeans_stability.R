#########################################################################################
# R script to run flowMeans for stability analysis
# 
# This script runs one iteration of flowMeans for the stability analysis. The main 
# script "stability_analysis.R" then runs the scripts for each clustering method several
# times in a loop.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the stability analysis.
# 
# Lukas M. Weber, February 2016
#########################################################################################


library(flowCore)
library(flowMeans)

# helper functions
source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for stability analysis
source("load_results_truth_stability.R")



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")

data_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_Levine_32, transformation = FALSE))
data_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_Levine_13, transformation = FALSE))
data_Nilsson <- flowCore::exprs(flowCore::read.FCS(file_Nilsson, transformation = FALSE))
data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

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


# subset data

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
data_Nilsson <- data_Nilsson[, marker_cols_Nilsson]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

# dim(data_Levine_32)
# dim(data_Levine_13)
# dim(data_Nilsson)
# dim(data_Mosmann)




#####################
### Run flowMeans ###
#####################

# run flowMeans

# set.seed(123)
runtime_Levine_32 <- system.time({
  out_flowMeans_Levine_32 <- flowMeans(data_Levine_32, Standardize = FALSE)
})


# set.seed(123)
runtime_Levine_13 <- system.time({
  out_flowMeans_Levine_13 <- flowMeans(data_Levine_13, Standardize = FALSE)
})


# set.seed(123)
runtime_Nilsson <- system.time({
  out_flowMeans_Nilsson <- flowMeans(data_Nilsson, Standardize = FALSE, NumC = 50)
})


# set.seed(123)
runtime_Mosmann <- system.time({
  out_flowMeans_Mosmann <- flowMeans(data_Mosmann, Standardize = FALSE, NumC = 50)
})


# extract cluster labels

clus_flowMeans_Levine_32 <- out_flowMeans_Levine_32@Label
clus_flowMeans_Levine_13 <- out_flowMeans_Levine_13@Label
clus_flowMeans_Nilsson <- out_flowMeans_Nilsson@Label
clus_flowMeans_Mosmann <- out_flowMeans_Mosmann@Label

# length(clus_flowMeans_Levine_32)
# length(clus_flowMeans_Levine_13)
# length(clus_flowMeans_Nilsson)
# length(clus_flowMeans_Mosmann)


# cluster sizes and number of clusters

# table(clus_flowMeans_Levine_32)
# table(clus_flowMeans_Levine_13)
# table(clus_flowMeans_Nilsson)
# table(clus_flowMeans_Mosmann)

# length(table(clus_flowMeans_Levine_32))
# length(table(clus_flowMeans_Levine_13))
# length(table(clus_flowMeans_Nilsson))
# length(table(clus_flowMeans_Mosmann))




#########################
### CALCULATE RESULTS ###
#########################

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_flowMeans_Levine_32 <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine_32, clus_truth_Levine_32)
res_flowMeans_Levine_13 <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine_13, clus_truth_Levine_13)
res_flowMeans_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Nilsson, clus_truth_Nilsson)
res_flowMeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Mosmann, clus_truth_Mosmann)

# res_flowMeans_Levine_32
# res_flowMeans_Levine_13
# res_flowMeans_Nilsson
# res_flowMeans_Mosmann




##############################################################
### OUTPUT RESULTS FOR ONE ITERATION OF STABILITY ANALYSIS ###
##############################################################

# data sets with multiple populations (Levine_32, Levine_13)
# output mean F1 score, mean precision, mean recall, runtime

res_stability_flowMeans_Levine_32 <- list(
  mean_F1 = mean(res_flowMeans_Levine_32$F1), 
  mean_pr = mean(res_flowMeans_Levine_32$pr), 
  mean_re = mean(res_flowMeans_Levine_32$re), 
  runtime = unname(runtime_Levine_32["elapsed"])
)

res_stability_flowMeans_Levine_13 <- list(
  mean_F1 = mean(res_flowMeans_Levine_13$F1), 
  mean_pr = mean(res_flowMeans_Levine_13$pr), 
  mean_re = mean(res_flowMeans_Levine_13$re), 
  runtime = unname(runtime_Levine_13["elapsed"])
)


# data sets with a single rare population of interest (Nilsson, Mosmann)
# output F1 score, precision, recall (for population of interest), and runtime

res_stability_flowMeans_Nilsson <- list(
  F1 = as.numeric(res_flowMeans_Nilsson$F1), 
  pr = as.numeric(res_flowMeans_Nilsson$pr), 
  re = as.numeric(res_flowMeans_Nilsson$re), 
  runtime = unname(runtime_Nilsson["elapsed"])
)

res_stability_flowMeans_Mosmann <- list(
  F1 = as.numeric(res_flowMeans_Mosmann$F1), 
  pr = as.numeric(res_flowMeans_Mosmann$pr), 
  re = as.numeric(res_flowMeans_Mosmann$re), 
  runtime = unname(runtime_Mosmann["elapsed"])
)


