#########################################################################################
# R script to run and evaluate FlowSOM for range of values for number of clusters k
#
# Lukas Weber, September 2016
#########################################################################################


### see scripts "run_FlowSOM.R" and "evaluate_FlowSOM.R" for more detailed code on how to
### run FlowSOM and calculate mean F1 scores


library(flowCore)
library(FlowSOM)




#################
### LOAD DATA ###
#################

# filenames

DATA_DIR <- "../../../benchmark_data_sets"

files <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik/data/Samusik_all.fcs"), 
  Nilsson_rare = file.path(DATA_DIR, "Nilsson_rare/data/Nilsson_rare.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare.fcs")
)

is_rare <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)


# load data files: FlowSOM requires flowFrame objects

data <- vector("list", length(files))
names(data) <- names(files)

for (i in 1:length(data)) {
  f <- files[[i]]
  data[[i]] <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
}

sapply(data, dim)


# indices of protein marker columns

marker_cols <- list(
  Levine_32dim = 5:36, 
  Levine_13dim = 1:13, 
  Samusik_01   = 9:47, 
  Samusik_all  = 9:47, 
  Nilsson_rare = c(5:7, 9:18), 
  Mosmann_rare = c(7:9, 11:21)
)

sapply(marker_cols, length)


# extract true population labels

clus_truth <- vector("list", length(files))
names(clus_truth) <- names(files)

for (i in 1:length(clus_truth)) {
  data_truth_i <- flowCore::exprs(flowCore::read.FCS(files[[i]], transformation = FALSE, truncate_max_range = FALSE))
  clus_truth[[i]] <- data_truth_i[, "label"]
}

sapply(clus_truth, length)
lapply(clus_truth, table)




#############################################################
### Run FlowSOM_pre: manually selected number of clusters ###
#############################################################

# run FlowSOM_pre (i.e. initial steps of FlowSOM algorithm) with manually selected number
# of clusters


# grid size 20x20 (400 clusters) for Mosmann_rare (data set with very rare population);
# and grid size 10x10 (100 clusters, i.e. default) for all other data sets

# grid sizes
grid_size <- list(
  Levine_32dim = 10, 
  Levine_13dim = 10, 
  Samusik_01   = 10, 
  Samusik_all  = 10, 
  Nilsson_rare = 10, 
  Mosmann_rare = 20
)

seed <- 100
out <- runtimes <- vector("list", length(data))
names(out) <- names(runtimes) <- names(data)

for (i in 1:length(data)) {
  set.seed(seed)
  fSOM <- FlowSOM::ReadInput(data[[i]], transform = FALSE, scale = FALSE)
  fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols[[i]], 
                            xdim = grid_size[[i]], ydim = grid_size[[i]])
  fSOM <- FlowSOM::BuildMST(fSOM)
  out[[i]] <- fSOM
}

# store output for meta-clustering step
out_pre_manual <- out




################################################################################
### Run FlowSOM (meta-clustering step): manually selected number of clusters ###
################################################################################

# run FlowSOM (meta-clustering step) for range of values of k, and calculate mean F1
# scores directly


library(clue)

# helper functions to match clusters and evaluate
source("../helpers/helper_match_evaluate_multiple.R")
source("../helpers/helper_match_evaluate_single.R")


# choose values of k
k_range <- seq(5, 80, by = 5)


# run FlowSOM and evaluate

seed <- 100

clus_range_k <- res_range_k <- vector("list", length(out_pre_manual))
names(clus_range_k) <- names(res_range_k) <- names(out_pre_manual)

for (i in 1:length(out_pre_manual)) {
  
  clus_range_k[[i]] <- res_range_k[[i]] <- vector("list", length(k_range))
  names(clus_range_k[[i]]) <- names(res_range_k[[i]]) <- k_range
  
  for (j in 1:length(k_range)) {
    
    # run meta-clustering
    set.seed(seed)
    meta <- suppressMessages(
      ConsensusClusterPlus::ConsensusClusterPlus(t(out_pre_manual[[i]]$map$codes), maxK = k_range[j])
    )
    meta <- meta[[k_range[j]]]$consensusClass
    
    # extract and store cluster labels
    clus <- meta[out_pre_manual[[i]]$map$mapping[, 1]]
    clus_range_k[[i]][[j]] <- clus
    
    # calculate and store results (mean F1 scores / F1 scores)
    if (!is_rare[i]) {
      res <- helper_match_evaluate_multiple(clus_range_k[[i]][[j]], clus_truth[[i]])
      res_range_k[[i]][[j]] <- res[["mean_F1"]]
      
    } else if (is_rare[i]) {
      res <- helper_match_evaluate_single(clus_range_k[[i]][[j]], clus_truth[[i]])
      res_range_k[[i]][[j]] <- res[["F1"]]
    }
  }
}


# collapse into vectors

df_res_range_k <- lapply(res_range_k, function(vals) sapply(vals, t))


# save results

for (i in 1:length(df_res_range_k)) {
  if (!is_rare[i]) {
    res_save <- data.frame(k = k_range, mean_F1 = df_res_range_k[[i]])
  } else if (is_rare[i]) {
    res_save <- data.frame(k = k_range, F1 = df_res_range_k[[i]])
  }
  
  filename <- paste0("../../results_range_k/results_FlowSOM_range_k_", names(df_res_range_k)[i], ".txt")
  
  write.table(res_save, file = filename, row.names = FALSE, quote = FALSE, sep = "\t")
}


# save session information

sink(file = "../../results_range_k/session_info_FlowSOM_range_k.txt")
sessionInfo()
sink()



