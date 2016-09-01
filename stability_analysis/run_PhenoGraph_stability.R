#########################################################################################
# Stability analysis (multiple random starts):
# Function to run and evaluate PhenoGraph once for each data set
#
# Lukas Weber, September 2016
#########################################################################################


# Note: PhenoGraph is available as a Python3 library, so this function runs external 
# Python3 scripts from the command line. This requires Python3 to be installed, along 
# with the numpy and phenograph libraries. Intermediate results are saved to external 
# files and loaded back into R, for compatibility with the rest of the stability analysis
# scripts.


run_PhenoGraph_stability <- function(data) {
  
  py_scripts_PhenoGraph <- list(
    Levine_32dim = "run_PhenoGraph_Levine_32dim.py", 
    Mosmann_rare = "run_PhenoGraph_Mosmann_rare.py"
  )
  
  out_files_PhenoGraph <- list(
    Levine_32dim = "../../results_stability_random_starts/PhenoGraph/python_out_Levine_32dim.txt", 
    Mosmann_rare = "../../results_stability_random_starts/PhenoGraph/python_out_Mosmann_rare.txt"
  )
  
  # run PhenoGraph once for each data set
  for (i in 1:length(data)) {
    cmd <- paste("python3", py_scripts_PhenoGraph[[i]])
    system(cmd)
  }
  
  # extract cluster labels
  clus <- vector("list", length(data))
  names(clus) <- names(data)
  
  for (i in 1:length(clus)) {
    clus[[i]] <- unname(unlist(read.table(out_files_PhenoGraph[[i]], header = FALSE, sep = "\t")))
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


