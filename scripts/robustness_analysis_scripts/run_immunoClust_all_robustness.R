#########################################################################################
# R script to run immunoClust_all for robustness analysis
# 
# This script runs one iteration of immunoClust_all for the robustness analysis. The main
# script "robustness_analysis.R" then runs the scripts for each clustering method several
# times in a loop.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the robustness analysis.
# 
# Lukas M. Weber, February 2016
#########################################################################################


# note installation from Bioconductor requires GNU Scientific Library

library(flowCore)
library(immunoClust)

# helper functions
source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for robustness analysis
source("load_results_truth_robustness.R")



#################
### LOAD DATA ###
#################

# use non-transformed data files, since immunoClust will transform automatically

DATA_DIR <- "../../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransform.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13_notransform.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC_notransform.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ_notransform.fcs")


# input data as flowFrame objects

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


# column names (parameters)

pars_Levine_32 <- colnames(data_Levine_32)[marker_cols_Levine_32]
pars_Levine_13 <- colnames(data_Levine_13)[marker_cols_Levine_13]
pars_Nilsson <- colnames(data_Nilsson)[marker_cols_Nilsson]
pars_Mosmann <- colnames(data_Mosmann)[marker_cols_Mosmann]

# length(pars_Levine_32)
# length(pars_Levine_13)
# length(pars_Nilsson)
# length(pars_Mosmann)




###########################################
### Run immunoClust and immunoClust_all ###
###########################################

# run immunoClust
# (note: decreasing the bias argument increases the number of clusters)

# set.seed(123)
# runtime_immunoClust_Levine_32 <- system.time({
#   out_immunoClust_Levine_32 <- immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32)
# })

# set.seed(123)
# runtime_immunoClust_Levine_13 <- system.time({
#   out_immunoClust_Levine_13 <- immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13)
# })

# set.seed(123)
# runtime_immunoClust_Nilsson <- system.time({
#   out_immunoClust_Nilsson <- immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, bias = 0.1)
# })

# set.seed(123)
# runtime_immunoClust_Mosmann <- system.time({
#   out_immunoClust_Mosmann <- immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, bias = 0.1)
# })


# run immmunoClust_all (with additional step to classify all cells)

# set.seed(123)
runtime_immunoClust_all_Levine_32 <- system.time({
  out_immunoClust_all_Levine_32 <- immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32, 
                                                             classify.all = TRUE)
})

# set.seed(123)
runtime_immunoClust_all_Levine_13 <- system.time({
  out_immunoClust_all_Levine_13 <- immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13, 
                                                             classify.all = TRUE)
})

# set.seed(123)
runtime_immunoClust_all_Nilsson <- system.time({
  out_immunoClust_all_Nilsson <- immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, 
                                                           bias = 0.1, classify.all = TRUE)
})

# set.seed(123)
runtime_immunoClust_all_Mosmann <- system.time({
  out_immunoClust_all_Mosmann <- immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, 
                                                           bias = 0.1, classify.all = TRUE)
})


# extract cluster labels

# immunoClust

# summary(out_immunoClust_Levine_32)  # number of clusters
# summary(out_immunoClust_Levine_13)
# summary(out_immunoClust_Nilsson)
# summary(out_immunoClust_Mosmann)

# clus_immunoClust_Levine_32 <- out_immunoClust_Levine_32@label
# clus_immunoClust_Levine_13 <- out_immunoClust_Levine_13@label
# clus_immunoClust_Nilsson <- out_immunoClust_Nilsson@label
# clus_immunoClust_Mosmann <- out_immunoClust_Mosmann@label


# immunoClust_all

# summary(out_immunoClust_all_Levine_32)
# summary(out_immunoClust_all_Levine_13)
# summary(out_immunoClust_all_Nilsson)
# summary(out_immunoClust_all_Mosmann)

clus_immunoClust_all_Levine_32 <- out_immunoClust_all_Levine_32@label
clus_immunoClust_all_Levine_13 <- out_immunoClust_all_Levine_13@label
clus_immunoClust_all_Nilsson <- out_immunoClust_all_Nilsson@label
clus_immunoClust_all_Mosmann <- out_immunoClust_all_Mosmann@label


# cluster sizes and number of clusters

# immunoClust

# table(clus_immunoClust_Levine_32)
# table(clus_immunoClust_Levine_13)
# table(clus_immunoClust_Nilsson)
# table(clus_immunoClust_Mosmann)

# length(table(clus_immunoClust_Levine_32))
# length(table(clus_immunoClust_Levine_13))
# length(table(clus_immunoClust_Nilsson))
# length(table(clus_immunoClust_Mosmann))

# immunoClust_all

# table(clus_immunoClust_all_Levine_32)
# table(clus_immunoClust_all_Levine_13)
# table(clus_immunoClust_all_Nilsson)
# table(clus_immunoClust_all_Mosmann)

# length(table(clus_immunoClust_all_Levine_32))
# length(table(clus_immunoClust_all_Levine_13))
# length(table(clus_immunoClust_all_Nilsson))
# length(table(clus_immunoClust_all_Mosmann))


# plots

# immunoClust

# data_transf_Levine_32 <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_32, data_Levine_32)
# data_transf_Levine_13 <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_13, data_Levine_13)
# data_transf_Nilsson <- immunoClust::trans.ApplyToData(out_immunoClust_Nilsson, data_Nilsson)
# data_transf_Mosmann <- immunoClust::trans.ApplyToData(out_immunoClust_Mosmann, data_Mosmann)

# png("../results/immunoClust/plot_immunoClust_Levine_2015_marrow_32.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_Levine_32, data_transf_Levine_32, N = 1000)
# dev.off()

# png("../results/immunoClust/plot_immunoClust_Levine_2015_marrow_13.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_Levine_13, data_transf_Levine_13, N = 1000)
# dev.off()

# png("../results/immunoClust/plot_immunoClust_Nilsson_2013_HSC.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_Nilsson, data_transf_Nilsson, N = 1000)
# dev.off()

# png("../results/immunoClust/plot_immunoClust_Mosmann_2014_activ.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_Mosmann, data_transf_Mosmann, N = 1000)
# dev.off()


# immunoClust_all

# data_transf_all_Levine_32 <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_32, data_Levine_32)
# data_transf_all_Levine_13 <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_13, data_Levine_13)
# data_transf_all_Nilsson <- immunoClust::trans.ApplyToData(out_immunoClust_all_Nilsson, data_Nilsson)
# data_transf_all_Mosmann <- immunoClust::trans.ApplyToData(out_immunoClust_all_Mosmann, data_Mosmann)

# png("../results/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_32.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_all_Levine_32, data_transf_all_Levine_32, N = 1000)
# dev.off()

# png("../results/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_13.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_all_Levine_13, data_transf_all_Levine_13, N = 1000)
# dev.off()

# png("../results/immunoClust_all/plot_immunoClust_all_Nilsson_2013_HSC.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_all_Nilsson, data_transf_all_Nilsson, N = 1000)
# dev.off()

# png("../results/immunoClust_all/plot_immunoClust_all_Mosmann_2014_activ.png", width = 1000, height = 1000)
# immunoClust::splom(out_immunoClust_all_Mosmann, data_transf_all_Mosmann, N = 1000)
# dev.off()




#########################
### CALCULATE RESULTS ###
#########################

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_immunoClust_all_Levine_32 <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine_32, clus_truth_Levine_32)
res_immunoClust_all_Levine_13 <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine_13, clus_truth_Levine_13)
res_immunoClust_all_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Nilsson, clus_truth_Nilsson)
res_immunoClust_all_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Mosmann, clus_truth_Mosmann)

# res_immunoClust_all_Levine_32
# res_immunoClust_all_Levine_13
# res_immunoClust_all_Nilsson
# res_immunoClust_all_Mosmann




###############################################################
### OUTPUT RESULTS FOR ONE ITERATION OF ROBUSTNESS ANALYSIS ###
###############################################################

# data sets with multiple populations (Levine_32, Levine_13)
# output mean F1 score, mean precision, mean recall, runtime

res_robust_immunoClust_all_Levine_32 <- list(
  mean_F1 = mean(res_immunoClust_all_Levine_32$F1), 
  mean_pr = mean(res_immunoClust_all_Levine_32$pr), 
  mean_re = mean(res_immunoClust_all_Levine_32$re), 
  runtime = unname(runtime_immunoClust_all_Levine_32["elapsed"])
)

res_robust_immunoClust_all_Levine_13 <- list(
  mean_F1 = mean(res_immunoClust_all_Levine_13$F1), 
  mean_pr = mean(res_immunoClust_all_Levine_13$pr), 
  mean_re = mean(res_immunoClust_all_Levine_13$re), 
  runtime = unname(runtime_immunoClust_all_Levine_13["elapsed"])
)


# data sets with a single rare population of interest (Nilsson, Mosmann)
# output F1 score, precision, recall (for population of interest), and runtime

res_robust_immunoClust_all_Nilsson <- list(
  F1 = as.numeric(res_immunoClust_all_Nilsson$F1), 
  pr = as.numeric(res_immunoClust_all_Nilsson$pr), 
  re = as.numeric(res_immunoClust_all_Nilsson$re), 
  runtime = unname(runtime_immunoClust_all_Nilsson["elapsed"])
)

res_robust_immunoClust_all_Mosmann <- list(
  F1 = as.numeric(res_immunoClust_all_Mosmann$F1), 
  pr = as.numeric(res_immunoClust_all_Mosmann$pr), 
  re = as.numeric(res_immunoClust_all_Mosmann$re), 
  runtime = unname(runtime_immunoClust_all_Mosmann["elapsed"])
)


