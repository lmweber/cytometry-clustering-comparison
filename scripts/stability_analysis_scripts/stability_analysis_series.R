#########################################################################################
# R script to run stability analysis for methods in series
#
# This script runs clustering methods multiple times in series, which is required for 
# methods that save results in files instead of returning R objects. The remaining 
# methods are run in parallel in "stability_analysis_parallel.R".
# 
# Note: Set random seeds for each loop, but not for individual iterations.
#
# Lukas M. Weber, March 2016
#########################################################################################


# number of times to run each method
n <- 30

# directory to save results
RESULTS_DIR <- "../../results_stability_analysis"




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

# no R session information since FLOCK is a command-line program




##################
### PhenoGraph ###
##################

### run on Linux server


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

# no R session information since using Python command-line implementation of PhenoGraph



