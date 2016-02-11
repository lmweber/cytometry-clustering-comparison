#########################################################################################
# R script for robustness analysis
#
# This script runs each clustering method several times in a loop, returning 
# distributions of the main evaluation criteria (F1 score etc).
# 
# Note: Set a single random seed for each loop/method, but not for individual iterations.
#
# Lukas M. Weber, February 2016
#########################################################################################


# number of times to run each method
n <- 30

# directory to save results
RESULTS_DIR <- "../../results/robustness_analysis"




#############
### FLOCK ###
#############

# helper functions
source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for robustness analysis
source("load_results_truth_robustness.R")


# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))


# load results for FLOCK (see script "run_FLOCK_robustness.R"; FLOCK is not available as an R package)

FLOCK_DIR <- "../../../FLOCK"

for (i in 1:n) {
  
  # note: need to run this code here instead of in a separate script since FLOCK runs
  # externally from the command line
  
  # filenames
  file_FLOCK_Levine_32 <- paste0(FLOCK_DIR, "/robustness_analysis/Levine_2015_marrow_32/flock_results_Levine_2015_marrow_32_iteration_", i, ".txt")
  file_FLOCK_Levine_13 <- paste0(FLOCK_DIR, "/robustness_analysis/Levine_2015_marrow_13/flock_results_Levine_2015_marrow_13_iteration_", i, ".txt")
  file_FLOCK_Nilsson <- paste0(FLOCK_DIR, "/robustness_analysis/Nilsson_2013_HSC/flock_results_Nilsson_2013_HSC_iteration_", i, ".txt")
  file_FLOCK_Mosmann <- paste0(FLOCK_DIR, "/robustness_analysis/Mosmann_2014_activ/flock_results_Mosmann_2014_activ_iteration_", i, ".txt")
  
  file_runtime_FLOCK_Levine_32 <- paste0(FLOCK_DIR, "/robustness_analysis/Levine_2015_marrow_32/runtime_Levine_2015_marrow_32_iteration_", i, ".txt")
  file_runtime_FLOCK_Levine_13 <- paste0(FLOCK_DIR, "/robustness_analysis/Levine_2015_marrow_13/runtime_Levine_2015_marrow_13_iteration_", i, ".txt")
  file_runtime_FLOCK_Nilsson <- paste0(FLOCK_DIR, "/robustness_analysis/Nilsson_2013_HSC/runtime_Nilsson_2013_HSC_iteration_", i, ".txt")
  file_runtime_FLOCK_Mosmann <- paste0(FLOCK_DIR, "/robustness_analysis/Mosmann_2014_activ/runtime_Mosmann_2014_activ_iteration_", i, ".txt")
  
  # data
  data_FLOCK_Levine_32 <- read.table(file_FLOCK_Levine_32, header = TRUE, sep = "\t")
  data_FLOCK_Levine_13 <- read.table(file_FLOCK_Levine_13, header = TRUE, sep = "\t")
  data_FLOCK_Nilsson <- read.table(file_FLOCK_Nilsson, header = TRUE, sep = "\t")
  data_FLOCK_Mosmann <- read.table(file_FLOCK_Mosmann, header = TRUE, sep = "\t")
  
  clus_FLOCK_Levine_32 <- data_FLOCK_Levine_32[, "Population"]
  clus_FLOCK_Levine_13 <- data_FLOCK_Levine_13[, "Population"]
  clus_FLOCK_Nilsson <- data_FLOCK_Nilsson[, "Population"]
  clus_FLOCK_Mosmann <- data_FLOCK_Mosmann[, "Population"]
  
  # runtime
  runtime_FLOCK_Levine_32 <- read.table(file_runtime_FLOCK_Levine_32, header = TRUE, sep = "\t")
  runtime_FLOCK_Levine_13 <- read.table(file_runtime_FLOCK_Levine_13, header = TRUE, sep = "\t")
  runtime_FLOCK_Nilsson <- read.table(file_runtime_FLOCK_Nilsson, header = TRUE, sep = "\t")
  runtime_FLOCK_Mosmann <- read.table(file_runtime_FLOCK_Mosmann, header = TRUE, sep = "\t")
  
  # calculate results
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  res_FLOCK_Levine_32 <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine_32, clus_truth_Levine_32)
  res_FLOCK_Levine_13 <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine_13, clus_truth_Levine_13)
  res_FLOCK_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Nilsson, clus_truth_Nilsson)
  res_FLOCK_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Mosmann, clus_truth_Mosmann)
  
  # output results
  
  # data sets with multiple populations (Levine_32, Levine_13)
  # output mean F1 score, mean precision, mean recall, runtime
  res_robust_FLOCK_Levine_32 <- list(
    mean_F1 = mean(res_FLOCK_Levine_32$F1), 
    mean_pr = mean(res_FLOCK_Levine_32$pr), 
    mean_re = mean(res_FLOCK_Levine_32$re), 
    runtime = runtime_FLOCK_Levine_32[, "x"]
  )
  
  res_robust_FLOCK_Levine_13 <- list(
    mean_F1 = mean(res_FLOCK_Levine_13$F1), 
    mean_pr = mean(res_FLOCK_Levine_13$pr), 
    mean_re = mean(res_FLOCK_Levine_13$re), 
    runtime = runtime_FLOCK_Levine_13[, "x"]
  )
  
  # data sets with a single rare population of interest (Nilsson, Mosmann)
  # output F1 score, precision, recall (for population of interest), and runtime
  res_robust_FLOCK_Nilsson <- list(
    F1 = as.numeric(res_FLOCK_Nilsson$F1), 
    pr = as.numeric(res_FLOCK_Nilsson$pr), 
    re = as.numeric(res_FLOCK_Nilsson$re), 
    runtime = runtime_FLOCK_Nilsson[, "x"]
  )
  
  res_robust_FLOCK_Mosmann <- list(
    F1 = as.numeric(res_FLOCK_Mosmann$F1), 
    pr = as.numeric(res_FLOCK_Mosmann$pr), 
    re = as.numeric(res_FLOCK_Mosmann$re), 
    runtime = runtime_FLOCK_Mosmann[, "x"]
  )
  
  # store results
  res_robust_Levine_32[i, ] <- unlist(res_robust_FLOCK_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_FLOCK_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_FLOCK_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_FLOCK_Mosmann)
}


# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "FLOCK", "robustness_FLOCK_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "FLOCK", "robustness_FLOCK_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "FLOCK", "robustness_FLOCK_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "FLOCK", "robustness_FLOCK_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "FLOCK", "session_info_robustness_FLOCK.txt"))
sessionInfo()
sink()




#################
### flowMeans ###
#################

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_flowMeans_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_flowMeans_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_flowMeans_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_flowMeans_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_flowMeans_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "flowMeans", "robustness_flowMeans_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "flowMeans", "robustness_flowMeans_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "flowMeans", "robustness_flowMeans_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "flowMeans", "robustness_flowMeans_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "flowMeans", "session_info_robustness_flowMeans.txt"))
sessionInfo()
sink()




###############
### FlowSOM ###
###############

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_FlowSOM_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_FlowSOM_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_FlowSOM_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_FlowSOM_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_FlowSOM_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "FlowSOM", "robustness_FlowSOM_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "FlowSOM", "robustness_FlowSOM_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "FlowSOM", "robustness_FlowSOM_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "FlowSOM", "robustness_FlowSOM_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "FlowSOM", "session_info_robustness_FlowSOM.txt"))
sessionInfo()
sink()




####################
### FlowSOM_meta ###
####################

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_FlowSOM_meta_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_FlowSOM_meta_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_FlowSOM_meta_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_FlowSOM_meta_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_FlowSOM_meta_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "FlowSOM_meta", "robustness_FlowSOM_meta_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "FlowSOM_meta", "robustness_FlowSOM_meta_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "FlowSOM_meta", "robustness_FlowSOM_meta_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "FlowSOM_meta", "robustness_FlowSOM_meta_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "FlowSOM_meta", "session_info_robustness_FlowSOM_meta.txt"))
sessionInfo()
sink()




###################
### immunoClust ###
###################

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_immunoClust_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_immunoClust_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_immunoClust_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_immunoClust_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_immunoClust_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "immunoClust", "robustness_immunoClust_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "immunoClust", "robustness_immunoClust_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "immunoClust", "robustness_immunoClust_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "immunoClust", "robustness_immunoClust_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "immunoClust", "session_info_robustness_immunoClust.txt"))
sessionInfo()
sink()




#######################
### immunoClust_all ###
#######################

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_immunoClust_all_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_immunoClust_all_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_immunoClust_all_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_immunoClust_all_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_immunoClust_all_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "immunoClust_all", "robustness_immunoClust_all_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "immunoClust_all", "robustness_immunoClust_all_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "immunoClust_all", "robustness_immunoClust_all_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "immunoClust_all", "robustness_immunoClust_all_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "immunoClust_all", "session_info_robustness_immunoClust_all.txt"))
sessionInfo()
sink()




###################
### PhenoGraph ###
###################

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_PhenoGraph_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_PhenoGraph_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_PhenoGraph_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_PhenoGraph_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_PhenoGraph_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "PhenoGraph", "robustness_PhenoGraph_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "PhenoGraph", "robustness_PhenoGraph_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "PhenoGraph", "robustness_PhenoGraph_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "PhenoGraph", "robustness_PhenoGraph_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information (not necessary here since PhenoGraph is a Python package)

# sink(file = file.path(RESULTS_DIR, "PhenoGraph", "session_info_robustness_PhenoGraph.txt"))
# sessionInfo()
# sink()




###################
### SamSPECTRAL ###
###################

set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_SamSPECTRAL_robustness.R")
  
  res_robust_Levine_32[i, ] <- unlist(res_robust_SamSPECTRAL_Levine_32)
  res_robust_Levine_13[i, ] <- unlist(res_robust_SamSPECTRAL_Levine_13)
  res_robust_Nilsson[i, ] <- unlist(res_robust_SamSPECTRAL_Nilsson)
  res_robust_Mosmann[i, ] <- unlist(res_robust_SamSPECTRAL_Mosmann)
}

# save results to files

write.table(res_robust_Levine_32, file = file.path(RESULTS_DIR, "SamSPECTRAL", "robustness_SamSPECTRAL_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Levine_13, file = file.path(RESULTS_DIR, "SamSPECTRAL", "robustness_SamSPECTRAL_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Nilsson, file = file.path(RESULTS_DIR, "SamSPECTRAL", "robustness_SamSPECTRAL_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_robust_Mosmann, file = file.path(RESULTS_DIR, "SamSPECTRAL", "robustness_SamSPECTRAL_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "SamSPECTRAL", "session_info_robustness_SamSPECTRAL.txt"))
sessionInfo()
sink()



