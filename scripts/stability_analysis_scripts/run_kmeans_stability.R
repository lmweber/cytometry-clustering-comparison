#########################################################################################
# Function to run k-means for stability analysis
# 
# This script contains a function to run k-means once for each data set for the
# stability analysis. The main script "stability_analysis.R" then calls this function
# multiple times in parallel with BiocParallel::bplapply.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the stability analysis.
# 
# Lukas M. Weber, March 2016
#########################################################################################


run_kmeans_stability <- function(data_list) {
  
  # extract data_list
  data_Levine_32 <- data_list[[1]]
  data_Levine_13 <- data_list[[2]]
  data_Nilsson <- data_list[[3]]
  data_Mosmann <- data_list[[4]]
  
  # run k-means once for each data set (without random seeds)
  
  # number of clusters
  k_Levine_32 <- 40
  k_Levine_13 <- 40
  k_Nilsson <- 40
  k_Mosmann <- 40
  
  # run k-means
  # note: additional iterations required
  
  runtime_Levine_32 <- system.time({
    out_kmeans_Levine_32 <- kmeans(data_Levine_32, k_Levine_32, iter.max = 50)
  })
  
  runtime_Levine_13 <- system.time({
    out_kmeans_Levine_13 <- kmeans(data_Levine_13, k_Levine_13, iter.max = 50)
  })
  
  runtime_Nilsson <- system.time({
    out_kmeans_Nilsson <- kmeans(data_Nilsson, k_Nilsson, iter.max = 50)
  })
  
  runtime_Mosmann <- system.time({
    out_kmeans_Mosmann <- kmeans(data_Mosmann, k_Mosmann, iter.max = 50)
  })
  
  # extract cluster labels
  
  clus_kmeans_Levine_32 <- out_kmeans_Levine_32$cluster
  clus_kmeans_Levine_13 <- out_kmeans_Levine_13$cluster
  clus_kmeans_Nilsson <- out_kmeans_Nilsson$cluster
  clus_kmeans_Mosmann <- out_kmeans_Mosmann$cluster
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  
  res_kmeans_Levine_32 <- helper_match_clusters_and_evaluate(clus_kmeans_Levine_32, clus_truth_Levine_32)
  res_kmeans_Levine_13 <- helper_match_clusters_and_evaluate(clus_kmeans_Levine_13, clus_truth_Levine_13)
  res_kmeans_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_kmeans_Nilsson, clus_truth_Nilsson)
  res_kmeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_kmeans_Mosmann, clus_truth_Mosmann)
  
  # return results for one iteration of stability analysis
  
  # data sets with multiple populations (Levine_32, Levine_13)
  # return mean F1 score, mean precision, mean recall, runtime
  
  res_stability_kmeans_Levine_32 <- c(
    mean_F1 = mean(res_kmeans_Levine_32$F1), 
    mean_pr = mean(res_kmeans_Levine_32$pr), 
    mean_re = mean(res_kmeans_Levine_32$re), 
    runtime = unname(runtime_Levine_32["elapsed"])
  )
  
  res_stability_kmeans_Levine_13 <- c(
    mean_F1 = mean(res_kmeans_Levine_13$F1), 
    mean_pr = mean(res_kmeans_Levine_13$pr), 
    mean_re = mean(res_kmeans_Levine_13$re), 
    runtime = unname(runtime_Levine_13["elapsed"])
  )
  
  # data sets with a single rare population of interest (Nilsson, Mosmann)
  # return F1 score, precision, recall, runtime
  
  res_stability_kmeans_Nilsson <- c(
    F1 = as.numeric(res_kmeans_Nilsson$F1), 
    pr = as.numeric(res_kmeans_Nilsson$pr), 
    re = as.numeric(res_kmeans_Nilsson$re), 
    runtime = unname(runtime_Nilsson["elapsed"])
  )
  
  res_stability_kmeans_Mosmann <- c(
    F1 = as.numeric(res_kmeans_Mosmann$F1), 
    pr = as.numeric(res_kmeans_Mosmann$pr), 
    re = as.numeric(res_kmeans_Mosmann$re), 
    runtime = unname(runtime_Mosmann["elapsed"])
  )
  
  # return as list
  
  list(res_stability_kmeans_Levine_32, 
       res_stability_kmeans_Levine_13, 
       res_stability_kmeans_Nilsson, 
       res_stability_kmeans_Mosmann)
}


