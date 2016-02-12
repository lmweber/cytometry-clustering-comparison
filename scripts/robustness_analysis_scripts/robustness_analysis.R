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

### note: run on Mac instead of Linux server, since FLOCK C code did not compile on the server


set.seed(1234)

# data sets with multiple populations
res_robust_Levine_32 <- res_robust_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_robust_Nilsson <- res_robust_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_FLOCK_robustness.R")
  
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

# save session information (not relevant here since FLOCK is a command-line program)

# sink(file = file.path(RESULTS_DIR, "FLOCK", "session_info_robustness_FLOCK.txt"))
# sessionInfo()
# sink()




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

# save session information (not relevant here since PhenoGraph is a Python package)

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



