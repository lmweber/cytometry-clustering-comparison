#########################################################################################
# R script for stability analysis
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
RESULTS_DIR <- "../../results/stability_analysis"




#############
### FLOCK ###
#############

### note: run on Mac instead of Linux server, since FLOCK C code did not compile on the server


set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_FLOCK_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_FLOCK_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_FLOCK_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_FLOCK_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_FLOCK_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "FLOCK", "stability_FLOCK_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "FLOCK", "stability_FLOCK_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "FLOCK", "stability_FLOCK_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "FLOCK", "stability_FLOCK_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information (not relevant here since FLOCK is a command-line program)

# sink(file = file.path(RESULTS_DIR, "FLOCK", "session_info_stability_FLOCK.txt"))
# sessionInfo()
# sink()




#################
### flowMeans ###
#################

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_flowMeans_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_flowMeans_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_flowMeans_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_flowMeans_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_flowMeans_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "flowMeans", "session_info_stability_flowMeans.txt"))
sessionInfo()
sink()




###############
### FlowSOM ###
###############

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_FlowSOM_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_FlowSOM_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_FlowSOM_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_FlowSOM_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_FlowSOM_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "FlowSOM", "session_info_stability_FlowSOM.txt"))
sessionInfo()
sink()




####################
### FlowSOM_meta ###
####################

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_FlowSOM_meta_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_FlowSOM_meta_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_FlowSOM_meta_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_FlowSOM_meta_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_FlowSOM_meta_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "FlowSOM_meta", "session_info_stability_FlowSOM_meta.txt"))
sessionInfo()
sink()




###################
### immunoClust ###
###################

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_immunoClust_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_immunoClust_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_immunoClust_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_immunoClust_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_immunoClust_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "immunoClust", "session_info_stability_immunoClust.txt"))
sessionInfo()
sink()




#######################
### immunoClust_all ###
#######################

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_immunoClust_all_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_immunoClust_all_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_immunoClust_all_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_immunoClust_all_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_immunoClust_all_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "immunoClust_all", "session_info_stability_immunoClust_all.txt"))
sessionInfo()
sink()




##################
### PhenoGraph ###
##################

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_PhenoGraph_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_PhenoGraph_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_PhenoGraph_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_PhenoGraph_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_PhenoGraph_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "PhenoGraph", "stability_PhenoGraph_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "PhenoGraph", "stability_PhenoGraph_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "PhenoGraph", "stability_PhenoGraph_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "PhenoGraph", "stability_PhenoGraph_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information (not relevant here since PhenoGraph is a Python package)

# sink(file = file.path(RESULTS_DIR, "PhenoGraph", "session_info_stability_PhenoGraph.txt"))
# sessionInfo()
# sink()




###################
### SamSPECTRAL ###
###################

set.seed(1234)

# data sets with multiple populations
res_stability_Levine_32 <- res_stability_Levine_13 <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("mean_F1", "mean_pr", "mean_re", "runtime")))

# data sets with one population of interest
res_stability_Nilsson <- res_stability_Mosmann <- 
  matrix(NA, nrow = n, ncol = 4, dimnames = list(NULL, c("F1", "pr", "re", "runtime")))

# run clustering method n times

for (i in 1:n) {
  source("run_SamSPECTRAL_stability.R")
  
  res_stability_Levine_32[i, ] <- unlist(res_stability_SamSPECTRAL_Levine_32)
  res_stability_Levine_13[i, ] <- unlist(res_stability_SamSPECTRAL_Levine_13)
  res_stability_Nilsson[i, ] <- unlist(res_stability_SamSPECTRAL_Nilsson)
  res_stability_Mosmann[i, ] <- unlist(res_stability_SamSPECTRAL_Mosmann)
}

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RESULTS_DIR, "SamSPECTRAL", "session_info_stability_SamSPECTRAL.txt"))
sessionInfo()
sink()



