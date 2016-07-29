#########################################################################################
# Function to run FlowSOM_meta for stability analysis
# 
# This script contains a function to run FlowSOM_meta once for each data set for the
# stability analysis. The main script "stability_analysis.R" then calls this function
# multiple times in parallel with BiocParallel::bplapply.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the stability analysis.
# 
# Lukas M. Weber, March 2016
#########################################################################################


run_FlowSOM_meta_stability <- function(data_list) {
  
  # extract data_list
  data_Levine_32 <- data_list[[1]]
  data_Levine_13 <- data_list[[2]]
  data_Nilsson <- data_list[[3]]
  data_Mosmann <- data_list[[4]]
  
  # run FlowSOM_meta once for each data set (without random seeds)
  
  # grid size (e.g. 10x10 or 20x20 grid, i.e. 100 or 400 clusters)
  grid_Levine_32 <- 10
  grid_Levine_13 <- 10
  grid_Nilsson <- 10
  grid_Mosmann <- 20
  
  # number of clusters
  k_Levine_32 <- 40
  k_Levine_13 <- 40
  k_Nilsson <- 40
  k_Mosmann <- 40
  
  # step 1: run FlowSOM
  
  runtime_FlowSOM_Levine_32 <- system.time({
    fSOM_Levine_32 <- FlowSOM::ReadInput(data_Levine_32, transform = FALSE, scale = FALSE)
    fSOM_Levine_32 <- FlowSOM::BuildSOM(fSOM_Levine_32, colsToUse = marker_cols_Levine_32, 
                                        xdim = grid_Levine_32, ydim = grid_Levine_32)
    fSOM_Levine_32 <- FlowSOM::BuildMST(fSOM_Levine_32)
  })
  
  runtime_FlowSOM_Levine_13 <- system.time({
    fSOM_Levine_13 <- FlowSOM::ReadInput(data_Levine_13, transform = FALSE, scale = FALSE)
    fSOM_Levine_13 <- FlowSOM::BuildSOM(fSOM_Levine_13, colsToUse = marker_cols_Levine_13, 
                                        xdim = grid_Levine_13, ydim = grid_Levine_13)
    fSOM_Levine_13 <- FlowSOM::BuildMST(fSOM_Levine_13)
  })
  
  runtime_FlowSOM_Nilsson <- system.time({
    fSOM_Nilsson <- FlowSOM::ReadInput(data_Nilsson, transform = FALSE, scale = FALSE)
    fSOM_Nilsson <- FlowSOM::BuildSOM(fSOM_Nilsson, colsToUse = marker_cols_Nilsson, 
                                      xdim = grid_Nilsson, ydim = grid_Nilsson)
    fSOM_Nilsson <- FlowSOM::BuildMST(fSOM_Nilsson)
  })
  
  runtime_FlowSOM_Mosmann <- system.time({
    fSOM_Mosmann <- FlowSOM::ReadInput(data_Mosmann, transform = FALSE, scale = FALSE)
    fSOM_Mosmann <- FlowSOM::BuildSOM(fSOM_Mosmann, colsToUse = marker_cols_Mosmann, 
                                      xdim = grid_Mosmann, ydim = grid_Mosmann)
    fSOM_Mosmann <- FlowSOM::BuildMST(fSOM_Mosmann)
  })
  
  # step 2: run meta-clustering
  
  runtime_FlowSOM_meta_Levine_32 <- system.time({
    meta_clustering_Levine_32 <- FlowSOM::metaClustering_consensus(fSOM_Levine_32$map$codes, k = k_Levine_32)
  })
  
  runtime_FlowSOM_meta_Levine_13 <- system.time({
    meta_clustering_Levine_13 <- FlowSOM::metaClustering_consensus(fSOM_Levine_13$map$codes, k = k_Levine_13)
  })
  
  runtime_FlowSOM_meta_Nilsson <- system.time({
    meta_clustering_Nilsson <- FlowSOM::metaClustering_consensus(fSOM_Nilsson$map$codes, k = k_Nilsson)
  })
  
  runtime_FlowSOM_meta_Mosmann <- system.time({
    meta_clustering_Mosmann <- FlowSOM::metaClustering_consensus(fSOM_Mosmann$map$codes, k = k_Mosmann)
  })
  
  # combine runtime
  
  runtime_FlowSOM_meta_Levine_32 <- runtime_FlowSOM_Levine_32 + runtime_FlowSOM_meta_Levine_32
  runtime_FlowSOM_meta_Levine_13 <- runtime_FlowSOM_Levine_13 + runtime_FlowSOM_meta_Levine_13
  runtime_FlowSOM_meta_Nilsson <- runtime_FlowSOM_Nilsson + runtime_FlowSOM_meta_Nilsson
  runtime_FlowSOM_meta_Mosmann <- runtime_FlowSOM_Mosmann + runtime_FlowSOM_meta_Mosmann
  
  # extract cluster labels
  
  clus_FlowSOM_meta_Levine_32 <- meta_clustering_Levine_32[fSOM_Levine_32$map$mapping[, 1]]
  clus_FlowSOM_meta_Levine_13 <- meta_clustering_Levine_13[fSOM_Levine_13$map$mapping[, 1]]
  clus_FlowSOM_meta_Nilsson <- meta_clustering_Nilsson[fSOM_Nilsson$map$mapping[, 1]]
  clus_FlowSOM_meta_Mosmann <- meta_clustering_Mosmann[fSOM_Mosmann$map$mapping[, 1]]
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  
  res_FlowSOM_meta_Levine_32 <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine_32, clus_truth_Levine_32)
  res_FlowSOM_meta_Levine_13 <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine_13, clus_truth_Levine_13)
  res_FlowSOM_meta_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Nilsson, clus_truth_Nilsson)
  res_FlowSOM_meta_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Mosmann, clus_truth_Mosmann)
  
  # return results for one iteration of stability analysis
  
  # data sets with multiple populations (Levine_32, Levine_13)
  # return mean F1 score, mean precision, mean recall, runtime
  
  res_stability_FlowSOM_meta_Levine_32 <- c(
    mean_F1 = mean(res_FlowSOM_meta_Levine_32$F1), 
    mean_pr = mean(res_FlowSOM_meta_Levine_32$pr), 
    mean_re = mean(res_FlowSOM_meta_Levine_32$re), 
    runtime = unname(runtime_FlowSOM_meta_Levine_32["elapsed"])
  )
  
  res_stability_FlowSOM_meta_Levine_13 <- c(
    mean_F1 = mean(res_FlowSOM_meta_Levine_13$F1), 
    mean_pr = mean(res_FlowSOM_meta_Levine_13$pr), 
    mean_re = mean(res_FlowSOM_meta_Levine_13$re), 
    runtime = unname(runtime_FlowSOM_meta_Levine_13["elapsed"])
  )
  
  # data sets with a single rare population of interest (Nilsson, Mosmann)
  # return F1 score, precision, recall, runtime
  
  res_stability_FlowSOM_meta_Nilsson <- c(
    F1 = as.numeric(res_FlowSOM_meta_Nilsson$F1), 
    pr = as.numeric(res_FlowSOM_meta_Nilsson$pr), 
    re = as.numeric(res_FlowSOM_meta_Nilsson$re), 
    runtime = unname(runtime_FlowSOM_meta_Nilsson["elapsed"])
  )
  
  res_stability_FlowSOM_meta_Mosmann <- c(
    F1 = as.numeric(res_FlowSOM_meta_Mosmann$F1), 
    pr = as.numeric(res_FlowSOM_meta_Mosmann$pr), 
    re = as.numeric(res_FlowSOM_meta_Mosmann$re), 
    runtime = unname(runtime_FlowSOM_meta_Mosmann["elapsed"])
  )
  
  # return as list
  
  list(res_stability_FlowSOM_meta_Levine_32, 
       res_stability_FlowSOM_meta_Levine_13, 
       res_stability_FlowSOM_meta_Nilsson, 
       res_stability_FlowSOM_meta_Mosmann)
}


