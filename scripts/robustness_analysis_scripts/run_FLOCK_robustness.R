#########################################################################################
# R script to run FLOCK for robustness analysis
# 
# This script runs one iteration of FLOCK for the robustness analysis. Note that FLOCK 
# runs from the command line (not available as an R package). The main script 
# "robustness_analysis.R" then runs the scripts for each clustering method several times 
# in a loop.
# 
# Lukas M. Weber, February 2016
#########################################################################################


### run this script from the FLOCK program directory
### FLOCK automatically saves results in file "flock_results.txt" in the current directory


source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for robustness analysis
source("load_results_truth_robustness.R")


ORIGINAL_DIR <- getwd()
FLOCK_DIR <- "../../../FLOCK"



##################
### DATA FILES ###
##################

# data files need to be subset already (marker columns only)

file_Levine_32 <- "Levine_2015_marrow_32_markers_only.txt"
file_Levine_13 <- "Levine_2015_marrow_13_markers_only.txt"
file_Nilsson <- "Nilsson_2013_HSC_markers_only.txt"
file_Mosmann <- "Mosmann_2014_activ_markers_only.txt"




#################
### RUN FLOCK ###
#################

# run FLOCK once for each data set and read in results as R objects (FLOCK overwrites the
# results file "flock_results.txt" each time)

setwd(FLOCK_DIR)


runtime_Levine_32 <- system.time(
  system(paste0("./flock2 ./", file_Levine_32))
)
res_FLOCK_Levine_32 <- read.table("flock_results.txt", header = TRUE, sep = "\t")
clus_FLOCK_Levine_32 <- res_FLOCK_Levine_32[, "Population"]


runtime_Levine_13 <- system.time(
  system(paste0("./flock2 ./", file_Levine_13))
)
res_FLOCK_Levine_13 <- read.table("flock_results.txt", header = TRUE, sep = "\t")
clus_FLOCK_Levine_13 <- res_FLOCK_Levine_13[, "Population"]


runtime_Nilsson <- system.time(
  system(paste0("./flock2 ./", file_Nilsson))
)
res_FLOCK_Nilsson <- read.table("flock_results.txt", header = TRUE, sep = "\t")
clus_FLOCK_Nilsson <- res_FLOCK_Nilsson[, "Population"]


runtime_Mosmann <- system.time(
  system(paste0("./flock2 ./", file_Mosmann))
)
res_FLOCK_Mosmann <- read.table("flock_results.txt", header = TRUE, sep = "\t")
clus_FLOCK_Mosmann <- res_FLOCK_Mosmann[, "Population"]


setwd(ORIGINAL_DIR)




#########################
### CALCULATE RESULTS ###
#########################

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FLOCK_Levine_32 <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine_32, clus_truth_Levine_32)
res_FLOCK_Levine_13 <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine_13, clus_truth_Levine_13)
res_FLOCK_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Nilsson, clus_truth_Nilsson)
res_FLOCK_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Mosmann, clus_truth_Mosmann)

# res_FLOCK_Levine_32
# res_FLOCK_Levine_13
# res_FLOCK_Nilsson
# res_FLOCK_Mosmann




###############################################################
### OUTPUT RESULTS FOR ONE ITERATION OF ROBUSTNESS ANALYSIS ###
###############################################################

# data sets with multiple populations (Levine_32, Levine_13)
# output mean F1 score, mean precision, mean recall, runtime

res_robust_FLOCK_Levine_32 <- list(
  mean_F1 = mean(res_FLOCK_Levine_32$F1), 
  mean_pr = mean(res_FLOCK_Levine_32$pr), 
  mean_re = mean(res_FLOCK_Levine_32$re), 
  runtime = unname(runtime_Levine_32["elapsed"])
)

res_robust_FLOCK_Levine_13 <- list(
  mean_F1 = mean(res_FLOCK_Levine_13$F1), 
  mean_pr = mean(res_FLOCK_Levine_13$pr), 
  mean_re = mean(res_FLOCK_Levine_13$re), 
  runtime = unname(runtime_Levine_13["elapsed"])
)


# data sets with a single rare population of interest (Nilsson, Mosmann)
# output F1 score, precision, recall (for population of interest), and runtime

res_robust_FLOCK_Nilsson <- list(
  F1 = as.numeric(res_FLOCK_Nilsson$F1), 
  pr = as.numeric(res_FLOCK_Nilsson$pr), 
  re = as.numeric(res_FLOCK_Nilsson$re), 
  runtime = unname(runtime_Nilsson["elapsed"])
)

res_robust_FLOCK_Mosmann <- list(
  F1 = as.numeric(res_FLOCK_Mosmann$F1), 
  pr = as.numeric(res_FLOCK_Mosmann$pr), 
  re = as.numeric(res_FLOCK_Mosmann$re), 
  runtime = unname(runtime_Mosmann["elapsed"])
)



