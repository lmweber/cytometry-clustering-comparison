#########################################################################################
# Function to run FlowSOM for stability analysis
# 
# This script contains a function to run FlowSOM once for each data set for the
# stability analysis. The main script "stability_analysis.R" then calls this function
# multiple times in parallel with BiocParallel::bplapply.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the stability analysis.
# 
# Lukas M. Weber, March 2016
#########################################################################################


run_FlowSOM_stability <- function(data_list) {
  
  # extract data_list
  data_Levine_32 <- data_list[[1]]
  data_Levine_13 <- data_list[[2]]
  data_Nilsson <- data_list[[3]]
  data_Mosmann <- data_list[[4]]
  
  # run FlowSOM once for each data set (without random seeds)
  
  # grid size (e.g. 10x10 or 20x20 grid, i.e. 100 or 400 clusters)
  grid_Levine_32 <- 10
  grid_Levine_13 <- 10
  grid_Nilsson <- 10
  grid_Mosmann <- 20
  
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
  
  # extract cluster labels
  
  clus_FlowSOM_Levine_32 <- fSOM_Levine_32$map$mapping[, 1]
  clus_FlowSOM_Levine_13 <- fSOM_Levine_13$map$mapping[, 1]
  clus_FlowSOM_Nilsson <- fSOM_Nilsson$map$mapping[, 1]
  clus_FlowSOM_Mosmann <- fSOM_Mosmann$map$mapping[, 1]
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  
  res_FlowSOM_Levine_32 <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine_32, clus_truth_Levine_32)
  res_FlowSOM_Levine_13 <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine_13, clus_truth_Levine_13)
  res_FlowSOM_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Nilsson, clus_truth_Nilsson)
  res_FlowSOM_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Mosmann, clus_truth_Mosmann)
  
  
  # return results for one iteration of stability analysis
  
  # data sets with multiple populations (Levine_32, Levine_13)
  # return mean F1 score, mean precision, mean recall, runtime
  
  res_stability_FlowSOM_Levine_32 <- c(
    mean_F1 = mean(res_FlowSOM_Levine_32$F1), 
    mean_pr = mean(res_FlowSOM_Levine_32$pr), 
    mean_re = mean(res_FlowSOM_Levine_32$re), 
    runtime = unname(runtime_FlowSOM_Levine_32["elapsed"])
  )
  
  res_stability_FlowSOM_Levine_13 <- c(
    mean_F1 = mean(res_FlowSOM_Levine_13$F1), 
    mean_pr = mean(res_FlowSOM_Levine_13$pr), 
    mean_re = mean(res_FlowSOM_Levine_13$re), 
    runtime = unname(runtime_FlowSOM_Levine_13["elapsed"])
  )
  
  # data sets with a single rare population of interest (Nilsson, Mosmann)
  # return F1 score, precision, recall, runtime
  
  res_stability_FlowSOM_Nilsson <- c(
    F1 = as.numeric(res_FlowSOM_Nilsson$F1), 
    pr = as.numeric(res_FlowSOM_Nilsson$pr), 
    re = as.numeric(res_FlowSOM_Nilsson$re), 
    runtime = unname(runtime_FlowSOM_Nilsson["elapsed"])
  )
  
  res_stability_FlowSOM_Mosmann <- c(
    F1 = as.numeric(res_FlowSOM_Mosmann$F1), 
    pr = as.numeric(res_FlowSOM_Mosmann$pr), 
    re = as.numeric(res_FlowSOM_Mosmann$re), 
    runtime = unname(runtime_FlowSOM_Mosmann["elapsed"])
  )
  
  # return as list
  
  list(res_stability_FlowSOM_Levine_32, 
       res_stability_FlowSOM_Levine_13, 
       res_stability_FlowSOM_Nilsson, 
       res_stability_FlowSOM_Mosmann)
}


