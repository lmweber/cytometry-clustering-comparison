#########################################################################################
# Function to run flowMeans for stability analysis
# 
# This script contains a function to run flowMeans once for each data set for the
# stability analysis. The main script "stability_analysis.R" then calls this function
# multiple times in parallel with BiocParallel::bplapply.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the stability analysis.
# 
# Lukas M. Weber, March 2016
#########################################################################################


run_flowMeans_stability <- function(data_list) {
  
  # extract data_list
  data_Levine_32 <- data_list[[1]]
  data_Levine_13 <- data_list[[2]]
  data_Nilsson <- data_list[[3]]
  data_Mosmann <- data_list[[4]]
  
  # run flowMeans once for each data set (without random seeds)
  
  # number of clusters
  k_Levine_32 <- 40
  k_Levine_13 <- 40
  k_Nilsson <- 40
  k_Mosmann <- 40
  
  runtime_Levine_32 <- system.time({
    out_flowMeans_Levine_32 <- flowMeans(data_Levine_32, Standardize = FALSE, NumC = k_Levine_32)
  })
  
  runtime_Levine_13 <- system.time({
    out_flowMeans_Levine_13 <- flowMeans(data_Levine_13, Standardize = FALSE, NumC = k_Levine_13)
  })
  
  runtime_Nilsson <- system.time({
    out_flowMeans_Nilsson <- flowMeans(data_Nilsson, Standardize = FALSE, NumC = k_Nilsson)
  })
  
  runtime_Mosmann <- system.time({
    out_flowMeans_Mosmann <- flowMeans(data_Mosmann, Standardize = FALSE, NumC = k_Mosmann)
  })
  
  # extract cluster labels
  
  clus_flowMeans_Levine_32 <- out_flowMeans_Levine_32@Label
  clus_flowMeans_Levine_13 <- out_flowMeans_Levine_13@Label
  clus_flowMeans_Nilsson <- out_flowMeans_Nilsson@Label
  clus_flowMeans_Mosmann <- out_flowMeans_Mosmann@Label
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  
  res_flowMeans_Levine_32 <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine_32, clus_truth_Levine_32)
  res_flowMeans_Levine_13 <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine_13, clus_truth_Levine_13)
  res_flowMeans_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Nilsson, clus_truth_Nilsson)
  res_flowMeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Mosmann, clus_truth_Mosmann)
  
  
  # return results for one iteration of stability analysis
  
  # data sets with multiple populations (Levine_32, Levine_13)
  # return mean F1 score, mean precision, mean recall, runtime
  
  res_stability_flowMeans_Levine_32 <- c(
    mean_F1 = mean(res_flowMeans_Levine_32$F1), 
    mean_pr = mean(res_flowMeans_Levine_32$pr), 
    mean_re = mean(res_flowMeans_Levine_32$re), 
    runtime = unname(runtime_Levine_32["elapsed"])
  )
  
  res_stability_flowMeans_Levine_13 <- c(
    mean_F1 = mean(res_flowMeans_Levine_13$F1), 
    mean_pr = mean(res_flowMeans_Levine_13$pr), 
    mean_re = mean(res_flowMeans_Levine_13$re), 
    runtime = unname(runtime_Levine_13["elapsed"])
  )
  
  # data sets with a single rare population of interest (Nilsson, Mosmann)
  # return F1 score, precision, recall, runtime
  
  res_stability_flowMeans_Nilsson <- c(
    F1 = as.numeric(res_flowMeans_Nilsson$F1), 
    pr = as.numeric(res_flowMeans_Nilsson$pr), 
    re = as.numeric(res_flowMeans_Nilsson$re), 
    runtime = unname(runtime_Nilsson["elapsed"])
  )
  
  res_stability_flowMeans_Mosmann <- c(
    F1 = as.numeric(res_flowMeans_Mosmann$F1), 
    pr = as.numeric(res_flowMeans_Mosmann$pr), 
    re = as.numeric(res_flowMeans_Mosmann$re), 
    runtime = unname(runtime_Mosmann["elapsed"])
  )
  
  # return as list
  
  list(res_stability_flowMeans_Levine_32, 
       res_stability_flowMeans_Levine_13, 
       res_stability_flowMeans_Nilsson, 
       res_stability_flowMeans_Mosmann)
}


