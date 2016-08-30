#########################################################################################
# Stability analysis (multiple random starts):
# Function to run and evaluate immunoClust once for each data set
#
# Lukas Weber, August 2016
#########################################################################################


random_starts_immmunoClust <- function(data) {
  
  # run once for each data set
  # note: don't set any random seeds, since we want a different random seed each time
  
  out <- vector("list", length(data))
  names(out) <- names(data)
  
  for (i in 1:length(out)) {
    data_i <- flowCore::flowFrame(data[[i]])  ## input data must be flowFrame
    out[[i]] <- immunoClust::cell.process(data[[i]], classify.all = TRUE)
  }
  
  # extract cluster labels
  clus <- vector("list", length(data))
  names(clus) <- names(data)
  
  for (i in 1:length(clus)) {
    clus[[i]] <- out[[i]]@label
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


