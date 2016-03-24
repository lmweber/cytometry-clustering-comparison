#########################################################################################
# R script to run FlowSOM_meta for range of values for number of clusters k
#
# Lukas M. Weber, March 2016
#########################################################################################


### see script "run_FlowSOM.R" for more detailed code to run FlowSOM / FlowSOM_meta


library(flowCore)
library(FlowSOM)



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")

# FlowSOM requires input data as flowFrame objects

data_Levine_32 <- flowCore::read.FCS(file_Levine_32, transformation = FALSE)
data_Levine_13 <- flowCore::read.FCS(file_Levine_13, transformation = FALSE)
data_Nilsson <- flowCore::read.FCS(file_Nilsson, transformation = FALSE)
data_Mosmann <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)

# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- c(7:9, 11:21)

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)




###################################
### Run initial steps (FlowSOM) ###
###################################

# run FlowSOM with manually selected number of clusters (10x10 or 20x20 grid)

# grid size (10x10 or 20x20 grid, i.e. 100 or 400 clusters)
grid_Levine_32 <- 10
grid_Levine_13 <- 10
grid_Nilsson <- 10
grid_Mosmann <- 20


set.seed(123)
runtime_FlowSOM_Levine_32 <- system.time({
  fSOM_Levine_32 <- FlowSOM::ReadInput(data_Levine_32, transform = FALSE, scale = FALSE)
  fSOM_Levine_32 <- FlowSOM::BuildSOM(fSOM_Levine_32, colsToUse = marker_cols_Levine_32, 
                                      xdim = grid_Levine_32, ydim = grid_Levine_32)
  fSOM_Levine_32 <- FlowSOM::BuildMST(fSOM_Levine_32)
})

set.seed(123)
runtime_FlowSOM_Levine_13 <- system.time({
  fSOM_Levine_13 <- FlowSOM::ReadInput(data_Levine_13, transform = FALSE, scale = FALSE)
  fSOM_Levine_13 <- FlowSOM::BuildSOM(fSOM_Levine_13, colsToUse = marker_cols_Levine_13, 
                                      xdim = grid_Levine_13, ydim = grid_Levine_13)
  fSOM_Levine_13 <- FlowSOM::BuildMST(fSOM_Levine_13)
})

set.seed(123)
runtime_FlowSOM_Nilsson <- system.time({
  fSOM_Nilsson <- FlowSOM::ReadInput(data_Nilsson, transform = FALSE, scale = FALSE)
  fSOM_Nilsson <- FlowSOM::BuildSOM(fSOM_Nilsson, colsToUse = marker_cols_Nilsson, 
                                    xdim = grid_Nilsson, ydim = grid_Nilsson)
  fSOM_Nilsson <- FlowSOM::BuildMST(fSOM_Nilsson)
})

set.seed(123)
runtime_FlowSOM_Mosmann <- system.time({
  fSOM_Mosmann <- FlowSOM::ReadInput(data_Mosmann, transform = FALSE, scale = FALSE)
  fSOM_Mosmann <- FlowSOM::BuildSOM(fSOM_Mosmann, colsToUse = marker_cols_Mosmann, 
                                    xdim = grid_Mosmann, ydim = grid_Mosmann)
  fSOM_Mosmann <- FlowSOM::BuildMST(fSOM_Mosmann)
})




##################################################################################
### Function to run metaclustering step (FlowSOM_meta) with variable number of ###
### clusters and calculate mean F1 score                                       ###
##################################################################################

run_FlowSOM_meta_k <- function(k, fSOM_Levine_32, fSOM_Levine_13, fSOM_Nilsson, fSOM_Mosmann) {
  
  # fSOM_Levine_32 etc are the outputs from the initial steps above (FlowSOM prior to metaclustering)
  
  # don't set random seeds inside function
  
  meta_clustering_Levine_32 <- FlowSOM::metaClustering_consensus(fSOM_Levine_32$map$codes, k = k)
  meta_clustering_Levine_13 <- FlowSOM::metaClustering_consensus(fSOM_Levine_13$map$codes, k = k)
  meta_clustering_Nilsson <- FlowSOM::metaClustering_consensus(fSOM_Nilsson$map$codes, k = k)
  meta_clustering_Mosmann <- FlowSOM::metaClustering_consensus(fSOM_Mosmann$map$codes, k = k)
  
  # extract cluster labels
  
  clus_FlowSOM_meta_Levine_32 <- meta_clustering_Levine_32[fSOM_Levine_32$map$mapping[, 1]]
  clus_FlowSOM_meta_Levine_13 <- meta_clustering_Levine_13[fSOM_Levine_13$map$mapping[, 1]]
  clus_FlowSOM_meta_Nilsson <- meta_clustering_Nilsson[fSOM_Nilsson$map$mapping[, 1]]
  clus_FlowSOM_meta_Mosmann <- meta_clustering_Mosmann[fSOM_Mosmann$map$mapping[, 1]]
  
  # match cluster labels by highest F1 score and calculate results
  # precision, recall, F1 score, matched cluster labels, number of cells per matched cluster
  
  res_FlowSOM_meta_Levine_32 <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine_32, clus_truth_Levine_32)
  res_FlowSOM_meta_Levine_13 <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine_13, clus_truth_Levine_13)
  res_FlowSOM_meta_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Nilsson, clus_truth_Nilsson)
  res_FlowSOM_meta_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Mosmann, clus_truth_Mosmann)
  
  # results: mean F1 score for data sets with multiple populations (Levine_32, Levine_13)
  
  res_k_Levine_32 <- c(
    mean_F1 = mean(res_FlowSOM_meta_Levine_32$F1), 
    mean_pr = mean(res_FlowSOM_meta_Levine_32$pr), 
    mean_re = mean(res_FlowSOM_meta_Levine_32$re)
  )
  
  res_k_Levine_13 <- c(
    mean_F1 = mean(res_FlowSOM_meta_Levine_13$F1), 
    mean_pr = mean(res_FlowSOM_meta_Levine_13$pr), 
    mean_re = mean(res_FlowSOM_meta_Levine_13$re)
  )
  
  # results: F1 score for data sets with a single rare population of interest (Nilsson, Mosmann)
  
  res_k_Nilsson <- c(
    F1 = as.numeric(res_FlowSOM_meta_Nilsson$F1), 
    pr = as.numeric(res_FlowSOM_meta_Nilsson$pr), 
    re = as.numeric(res_FlowSOM_meta_Nilsson$re)
  )
  
  res_k_Mosmann <- c(
    F1 = as.numeric(res_FlowSOM_meta_Mosmann$F1), 
    pr = as.numeric(res_FlowSOM_meta_Mosmann$pr), 
    re = as.numeric(res_FlowSOM_meta_Mosmann$re)
  )
  
  # return as list
  
  list(res_k_Levine_32, 
       res_k_Levine_13, 
       res_k_Nilsson, 
       res_k_Mosmann)
}




#################################################
### Run FlowSOM_meta for range of values of k ###
#################################################

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")

RES_DIR <- "../results_range_k"


# choose values of k

k_range <- seq(5, 80, by = 5)

# run FlowSOM_meta for range of values of k

set.seed(1234)
res_range_k <- lapply(k_range, run_FlowSOM_meta_k, 
                      fSOM_Levine_32, fSOM_Levine_13, fSOM_Nilsson, fSOM_Mosmann)

# collapse results into one data frame for each data set

res_range_k_Levine_32 <- cbind(k = k_range, t(sapply(res_range_k, function(s) s[[1]])))
res_range_k_Levine_13 <- cbind(k = k_range, t(sapply(res_range_k, function(s) s[[2]])))
res_range_k_Nilsson <- cbind(k = k_range, t(sapply(res_range_k, function(s) s[[3]])))
res_range_k_Mosmann <- cbind(k = k_range, t(sapply(res_range_k, function(s) s[[4]])))

# save results to files

write.table(res_range_k_Levine_32, file = file.path(RES_DIR, "results_FlowSOM_meta_range_k_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_range_k_Levine_13, file = file.path(RES_DIR, "results_FlowSOM_meta_range_k_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_range_k_Nilsson, file = file.path(RES_DIR, "results_FlowSOM_meta_range_k_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_range_k_Mosmann, file = file.path(RES_DIR, "results_FlowSOM_meta_range_k_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save session information

sink(file = file.path(RES_DIR, "session_FlowSOM_meta_range_k.txt"))
sessionInfo()
sink()



