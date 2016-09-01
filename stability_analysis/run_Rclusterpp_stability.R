#########################################################################################
# Stability analysis (multiple random starts):
# Function to run and evaluate Rclusterpp once for each data set
#
# Lukas Weber, September 2016
#########################################################################################


run_Rclusterpp_stability <- function(data) {
  
  # parameters
  
  n_cores <- 8
  
  k <- list(
    Levine_32dim = 40
  )
  
  method <- list(
    Levine_32dim = "ward"
  )
  
  # run once for each data set
  # note: don't set any random seeds, since we want a different random seed each time
  
  out <- vector("list", length(data))
  names(out) <- names(data)
  
  for (i in 1:length(out)) {
    Rclusterpp.setThreads(n_cores)
    out[[i]] <- Rclusterpp.hclust(data[[i]], method = method[[i]])
  }
  
  # extract cluster labels
  clus <- vector("list", length(data))
  names(clus) <- names(data)
  
  for (i in 1:length(clus)) {
    # cut dendrogram at k
    clus[[i]] <- cutree(out[[i]], k = k[[i]])
  }
  
  # calculate mean F1 scores / F1 scores
  res <- vector("list", length(clus))
  names(res) <- names(clus)
  
  for (i in 1:length(clus)) {
    if (!is_rare[i]) {
      res[[i]] <- helper_match_evaluate_multiple(clus[[i]], clus_truth[[i]])
      
    } else if (is_rare[i]) {
      res[[i]] <- helper_match_evaluate_single(clus[[i]], clus_truth[[i]])
    }
  }
  
  return(res)
}


