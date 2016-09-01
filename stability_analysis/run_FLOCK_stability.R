#########################################################################################
# Stability analysis:
# Function to run and evaluate FLOCK once for each data set
#
# Lukas Weber, September 2016
#########################################################################################


# note: cannot parallelize FLOCK, since it needs to read and write from files in the
# FLOCK program directory


run_FLOCK_stability <- function(data) {
  
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
  
  # run from FLOCK program directory
  CURRENT_DIR <- getwd()
  setwd("../../../algorithms/FLOCK")
  
  out <- vector("list", length(data))
  names(out) <- names(data)
  
  # run once for each data set
  for (i in 1:length(data)) {
    
    # save external data file
    write.table(data[[i]], file = "FLOCK_data_file.txt", quote = FALSE, sep = "\t", row.names = FALSE)
    
    # run FLOCK
    cmd <- "./flock2 ./FLOCK_data_file.txt"
    system(cmd)
    
    # read results from external results file
    out[[i]] <- read.table("flock_results.txt", header = TRUE, sep = "\t")
  }
  
  # extract cluster labels
  clus <- vector("list", length(data))
  names(clus) <- names(data)
  
  for (i in 1:length(clus)) {
    clus[[i]] <- out[[i]][, "Population"]
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
  
  # reset working directory
  setwd(CURRENT_DIR)
  
  return(res)
}


