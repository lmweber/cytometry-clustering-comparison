#########################################################################################
# Stability analysis (multiple random starts):
# Function to run and evaluate k-means once for each data set
#
# Lukas Weber, August 2016
#########################################################################################


random_starts_kmeans <- function(data) {
  
  # parameters
  
  # number of clusters k
  k <- list(
    Levine_32dim = 40, 
    Mosmann_rare = 40
  )
  
  iter.max <- 50  ## some random seeds require more iterations
  
  # run once for each data set
  # note: don't set any random seeds, since we want a different random seed each time
  
  out <- vector("list", length(data))
  names(out) <- names(data)
  
  for (i in 1:length(out)) {
    out[[i]] <- kmeans(data[[i]], centers = k[[i]], iter.max = iter.max)
  }
  
  # extract cluster labels
  clus <- vector("list", length(data))
  names(clus) <- names(data)
  
  for (i in 1:length(clus)) {
    clus[[i]] <- out[[i]]$cluster
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


