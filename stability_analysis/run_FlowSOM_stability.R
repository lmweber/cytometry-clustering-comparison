#########################################################################################
# Stability analysis (multiple random starts):
# Function to run and evaluate FlowSOM once for each data set
#
# Lukas Weber, September 2016
#########################################################################################


run_FlowSOM_stability <- function(data) {
  
  # parameters
  
  # grid sizes
  grid_size <- list(
    Levine_32dim = 10, 
    Mosmann_rare = 20
  )
  
  # number of clusters k
  k <- list(
    Levine_32dim = 40, 
    Mosmann_rare = 40
  )
  
  # run once for each data set
  # note: don't set any random seeds, since we want a different random seed each time
  
  out_pre <- out <- vector("list", length(data))
  names(out_pre) <- names(out) <- names(data)
  
  for (i in 1:length(out)) {
    data_i <- flowCore::flowFrame(data[[i]])  ## input data must be flowFrame
    fSOM <- FlowSOM::ReadInput(data_i, transform = FALSE, scale = FALSE)
    fSOM <- FlowSOM::BuildSOM(fSOM, 
                              colsToUse = NULL,  ## use all columns since already subsetted
                              xdim = grid_size[[i]], 
                              ydim = grid_size[[i]])
    #fSOM <- FlowSOM::BuildMST(fSOM)  ## not required
    out_pre[[i]] <- fSOM  ## need to store intermediate result to get cluster labels
    
    # meta-clustering step
    meta <- suppressMessages(
      ConsensusClusterPlus::ConsensusClusterPlus(t(fSOM$map$codes), maxK = k[[i]])
    )
    meta <- meta[[k[[i]]]]$consensusClass
    
    out[[i]] <- meta
  }
  
  # extract cluster labels
  clus <- vector("list", length(data))
  names(clus) <- names(data)
  
  for (i in 1:length(clus)) {
    clus[[i]] <- out[[i]][out_pre[[i]]$map$mapping[, 1]]
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


