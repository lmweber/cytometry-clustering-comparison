#########################################################################################
# Python3 script to run PhenoGraph for stability analysis
# 
# This script runs one iteration of PhenoGraph by calling the Python3 scripts 
# "run_PhenoGraph_python_Levine_32.py" etc (one for each data set) and loading results
# back into R. The main script "stability_analysis.R" then runs the scripts for each
# clustering method several times in a loop.
# 
# Lukas M. Weber, March 2016
#########################################################################################


### run this script from the stability analysis directory


source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for stability analysis
source("load_results_truth_stability.R")



######################
### RUN PHENOGRAPH ###
######################

# run external Python script to run one iteration of Phenograph and save results to text files

runtime_Levine_32 <- system.time(
  system("python3 run_PhenoGraph_python_Levine_32.py")
)

runtime_Levine_13 <- system.time(
  system("python3 run_PhenoGraph_python_Levine_13.py")
)

runtime_Nilsson <- system.time(
  system("python3 run_PhenoGraph_python_Nilsson.py")
)

runtime_Mosmann <- system.time(
  system("python3 run_PhenoGraph_python_Mosmann.py")
)




#################################
### LOAD PYTHON OUTPUT INTO R ###
#################################

# load cluster labels

OUT_DIR <- "../../results_stability_analysis/PhenoGraph"

file_out_Levine_32 <- file.path(OUT_DIR, "python_out_Levine_32.txt")
file_out_Levine_13 <- file.path(OUT_DIR, "python_out_Levine_13.txt")
file_out_Nilsson <- file.path(OUT_DIR, "python_out_Nilsson.txt")
file_out_Mosmann <- file.path(OUT_DIR, "python_out_Mosmann.txt")


# load cluster labels

clus_PhenoGraph_Levine_32 <- unname(unlist(read.table(file_out_Levine_32, header = FALSE, sep = "\t")))
clus_PhenoGraph_Levine_13 <- unname(unlist(read.table(file_out_Levine_13, header = FALSE, sep = "\t")))
clus_PhenoGraph_Nilsson <- unname(unlist(read.table(file_out_Nilsson, header = FALSE, sep = "\t")))
clus_PhenoGraph_Mosmann <- unname(unlist(read.table(file_out_Mosmann, header = FALSE, sep = "\t")))

# length(clus_PhenoGraph_Levine_32)
# length(clus_PhenoGraph_Levine_13)
# length(clus_PhenoGraph_Nilsson)
# length(clus_PhenoGraph_Mosmann)

# table(clus_PhenoGraph_Levine_32)
# table(clus_PhenoGraph_Levine_13)
# table(clus_PhenoGraph_Nilsson)
# table(clus_PhenoGraph_Mosmann)




#########################
### CALCULATE RESULTS ###
#########################

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_PhenoGraph_Levine_32 <- helper_match_clusters_and_evaluate(clus_PhenoGraph_Levine_32, clus_truth_Levine_32)
res_PhenoGraph_Levine_13 <- helper_match_clusters_and_evaluate(clus_PhenoGraph_Levine_13, clus_truth_Levine_13)
res_PhenoGraph_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_PhenoGraph_Nilsson, clus_truth_Nilsson)
res_PhenoGraph_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_PhenoGraph_Mosmann, clus_truth_Mosmann)

# res_PhenoGraph_Levine_32
# res_PhenoGraph_Levine_13
# res_PhenoGraph_Nilsson
# res_PhenoGraph_Mosmann




##############################################################
### OUTPUT RESULTS FOR ONE ITERATION OF STABILITY ANALYSIS ###
##############################################################

# data sets with multiple populations (Levine_32, Levine_13)
# output mean F1 score, mean precision, mean recall, runtime

res_stability_PhenoGraph_Levine_32 <- list(
  mean_F1 = mean(res_PhenoGraph_Levine_32$F1), 
  mean_pr = mean(res_PhenoGraph_Levine_32$pr), 
  mean_re = mean(res_PhenoGraph_Levine_32$re), 
  runtime = unname(runtime_Levine_32["elapsed"])
)

res_stability_PhenoGraph_Levine_13 <- list(
  mean_F1 = mean(res_PhenoGraph_Levine_13$F1), 
  mean_pr = mean(res_PhenoGraph_Levine_13$pr), 
  mean_re = mean(res_PhenoGraph_Levine_13$re), 
  runtime = unname(runtime_Levine_13["elapsed"])
)


# data sets with a single rare population of interest (Nilsson, Mosmann)
# output F1 score, precision, recall (for population of interest), and runtime

res_stability_PhenoGraph_Nilsson <- list(
  F1 = as.numeric(res_PhenoGraph_Nilsson$F1), 
  pr = as.numeric(res_PhenoGraph_Nilsson$pr), 
  re = as.numeric(res_PhenoGraph_Nilsson$re), 
  runtime = unname(runtime_Nilsson["elapsed"])
)

res_stability_PhenoGraph_Mosmann <- list(
  F1 = as.numeric(res_PhenoGraph_Mosmann$F1), 
  pr = as.numeric(res_PhenoGraph_Mosmann$pr), 
  re = as.numeric(res_PhenoGraph_Mosmann$re), 
  runtime = unname(runtime_Mosmann["elapsed"])
)



