#########################################################################################
# Stability analysis (multiple random starts):
# Function to run and evaluate FLOCK once for each data set
#
# Lukas Weber, August 2016
#########################################################################################


# note: cannot parallelize FLOCK, since it needs to read and write from files in the
# FLOCK program directory


random_starts_FLOCK <- function(data) {
  
  # run from FLOCK program directory
  CURRENT_DIR <- getwd()
  setwd("../../../algorithms/FLOCK")
  
  out <- runtimes <- vector("list", length(data))
  names(out) <- names(runtimes) <- names(data)
  
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


