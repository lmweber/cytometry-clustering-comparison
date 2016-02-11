#########################################################################################
# R script to run FLOCK for robustness analysis
# 
# This script runs n iterations of FLOCK from the command line (FLOCK is not available as
# an R package). The results are then loaded in the main script "robustness_analysis.R".
# 
# Lukas M. Weber, February 2016
#########################################################################################


### run from FLOCK program directory
### see script "run_FLOCK.R" for instructions on how to run FLOCK

### note: FLOCK automatically saves results in file "flock_results.txt" in the current 
### directory, so need to copy this file to a different directory after each iteration


# number of times to run
n <- 30

# data files
files <- c("Levine_2015_marrow_32_markers_only.txt", 
           "Levine_2015_marrow_13_markers_only.txt", 
           "Nilsson_2013_HSC_markers_only.txt", 
           "Mosmann_2014_activ_markers_only.txt")

dirs <- c("Levine_2015_marrow_32", 
          "Levine_2015_marrow_13", 
          "Nilsson_2013_HSC", 
          "Mosmann_2014_activ")


# run FLOCK n times and copy/rename results file after each iteration

for (i in 1:n) {
  for (j in 1:length(files)) {
    
    # run FLOCK
    cmd_run <- paste0("./flock2 ./", files[j])
    runtime <- system.time({
      system(cmd_run)
    })
    
    # copy results file
    cmd_copy <- paste0("cp flock_results.txt robustness_analysis/", dirs[j], 
                       "/flock_results_", dirs[j], "_iteration_", i, ".txt")
    system(cmd_copy)
    
    # save runtime
    write.table(runtime["elapsed"], 
                file = paste0("robustness_analysis/", dirs[j], "/runtime_", dirs[j], "_iteration_", i, ".txt"), 
                row.names = FALSE, quote = FALSE, sep = "\t")
  }
}


