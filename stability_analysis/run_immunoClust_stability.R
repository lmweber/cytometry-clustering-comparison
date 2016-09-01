#########################################################################################
# Stability analysis:
# Function to run and evaluate immunoClust once for each data set
#
# Lukas Weber, September 2016
#########################################################################################


run_immunoClust_stability <- function(data) {
  
  # extract true population labels
  clus_truth <- vector("list", length(data))
  names(clus_truth) <- names(data)
  
  for (i in 1:length(clus_truth)) {
    clus_truth[[i]] <- data[[i]][, "label"]
  }
  
  # subset data: protein marker columns only
  marker_cols <- list(
    Levine_32dim = 5:36, 
    Mosmann_rare = c(7:9, 11:21)
  )
  
  for (i in 1:length(data)) {
    data[[i]] <- data[[i]][, marker_cols[[i]]]
  }
  
  # run once for each data set
  # note: don't set any random seeds, since we want a different random seed each time
  out <- vector("list", length(data))
  names(out) <- names(data)
  
  for (i in 1:length(out)) {
    data_i <- flowCore::flowFrame(data[[i]])  ## input data must be flowFrame
    out[[i]] <- immunoClust::cell.process(data_i, classify.all = TRUE)
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


