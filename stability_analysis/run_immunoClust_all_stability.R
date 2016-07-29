#########################################################################################
# Function to run immunoClust_all for stability analysis
# 
# This script contains a function to run immunoClust_all once for each data set for the
# stability analysis. The main script "stability_analysis.R" then calls this function
# multiple times in parallel with BiocParallel::bplapply.
# 
# Note: Don't set any random seeds here, since we want a different random start for each 
# iteration in the stability analysis.
# 
# Lukas M. Weber, March 2016
#########################################################################################


run_immunoClust_all_stability <- function(data_list) {
  
  # extract data_list
  data_Levine_32 <- data_list[[1]]
  data_Levine_13 <- data_list[[2]]
  data_Nilsson <- data_list[[3]]
  data_Mosmann <- data_list[[4]]
  
  # run immunoClust_all once for each data set (without random seeds)
  # i.e. immunoClust with additional step to classify all cells
  
  runtime_Levine_32 <- system.time({
    out_immunoClust_all_Levine_32 <- 
      immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32, classify.all = TRUE)
  })
  
  runtime_Levine_13 <- system.time({
    out_immunoClust_all_Levine_13 <- 
      immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13, classify.all = TRUE)
  })
  
  runtime_Nilsson <- system.time({
    out_immunoClust_all_Nilsson <- 
      immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, bias = 0.1, classify.all = TRUE)
  })
  
  runtime_Mosmann <- system.time({
    out_immunoClust_all_Mosmann <- 
      immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, bias = 0.1, classify.all = TRUE)
  })
  
  # extract cluster labels
  
  clus_immunoClust_all_Levine_32 <- out_immunoClust_all_Levine_32@label
  clus_immunoClust_all_Levine_13 <- out_immunoClust_all_Levine_13@label
  clus_immunoClust_all_Nilsson <- out_immunoClust_all_Nilsson@label
  clus_immunoClust_all_Mosmann <- out_immunoClust_all_Mosmann@label
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  
  res_immunoClust_all_Levine_32 <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine_32, clus_truth_Levine_32)
  res_immunoClust_all_Levine_13 <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine_13, clus_truth_Levine_13)
  res_immunoClust_all_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Nilsson, clus_truth_Nilsson)
  res_immunoClust_all_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Mosmann, clus_truth_Mosmann)
  
  # return results for one iteration of stability analysis
  
  # data sets with multiple populations (Levine_32, Levine_13)
  # return mean F1 score, mean precision, mean recall, runtime
  
  res_stability_immunoClust_all_Levine_32 <- c(
    mean_F1 = mean(res_immunoClust_all_Levine_32$F1), 
    mean_pr = mean(res_immunoClust_all_Levine_32$pr), 
    mean_re = mean(res_immunoClust_all_Levine_32$re), 
    runtime = unname(runtime_Levine_32["elapsed"])
  )
  
  res_stability_immunoClust_all_Levine_13 <- c(
    mean_F1 = mean(res_immunoClust_all_Levine_13$F1), 
    mean_pr = mean(res_immunoClust_all_Levine_13$pr), 
    mean_re = mean(res_immunoClust_all_Levine_13$re), 
    runtime = unname(runtime_Levine_13["elapsed"])
  )
  
  # data sets with a single rare population of interest (Nilsson, Mosmann)
  # return F1 score, precision, recall, runtime
  
  res_stability_immunoClust_all_Nilsson <- c(
    F1 = as.numeric(res_immunoClust_all_Nilsson$F1), 
    pr = as.numeric(res_immunoClust_all_Nilsson$pr), 
    re = as.numeric(res_immunoClust_all_Nilsson$re), 
    runtime = unname(runtime_Nilsson["elapsed"])
  )
  
  res_stability_immunoClust_all_Mosmann <- c(
    F1 = as.numeric(res_immunoClust_all_Mosmann$F1), 
    pr = as.numeric(res_immunoClust_all_Mosmann$pr), 
    re = as.numeric(res_immunoClust_all_Mosmann$re), 
    runtime = unname(runtime_Mosmann["elapsed"])
  )
  
  # return as list
  
  list(res_stability_immunoClust_all_Levine_32, 
       res_stability_immunoClust_all_Levine_13, 
       res_stability_immunoClust_all_Nilsson, 
       res_stability_immunoClust_all_Mosmann)
}


