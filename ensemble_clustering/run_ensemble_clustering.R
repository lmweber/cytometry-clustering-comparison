#########################################################################################
# R script to run ensemble clustering
#
# Lukas Weber, September 2016
#########################################################################################


library(clue)


# load clustering results (runtime: 15 mins)
source("../evaluate_results/evaluate_ClusterX.R")
source("../evaluate_results/evaluate_DensVM.R")
source("../evaluate_results/evaluate_FLOCK.R")
source("../evaluate_results/evaluate_flowClust.R")
source("../evaluate_results/evaluate_flowMeans.R")
source("../evaluate_results/evaluate_flowPeaks.R")
source("../evaluate_results/evaluate_FlowSOM.R")
source("../evaluate_results/evaluate_immunoClust.R")
source("../evaluate_results/evaluate_kmeans.R")
source("../evaluate_results/evaluate_PhenoGraph.R")
source("../evaluate_results/evaluate_Rclusterpp.R")
source("../evaluate_results/evaluate_SamSPECTRAL.R")
source("../evaluate_results/evaluate_SPADE.R")
source("../evaluate_results/evaluate_SWIFT.R")
source("../evaluate_results/evaluate_Xshift.R")




###########################
### ENSEMBLE CLUSTERING ###
###########################

# calculate ensemble clustering using all methods for each data set, excluding the 
# following:
# - methods that use subsampled data, since consensus clustering requires same data for
# each method (see parameters spreadsheet)
# - methods that remove outliers, since consensus clustering requires same data for each 
# method (SamSPECTRAL for all data sets; flowClust for Nilsson_rare and Mosmann_rare;
# SWIFT for data sets with multiple populations; X-shift for all data sets except
# Samusik_01 and Nilsson_rare)
# - methods that give a large number of small clusters (FlowSOM_pre; SWIFT for data sets
# with multiple populations), since this greatly slows down runtime


# create partitions in format required by CLUE

partition_Levine_32dim <- list(as.cl_partition(clus_FLOCK[[1]]), 
                               as.cl_partition(clus_flowMeans[[1]]), 
                               as.cl_partition(clus_flowPeaks[[1]]), 
                               as.cl_partition(clus_FlowSOM[[1]]), 
                               as.cl_partition(clus_kmeans[[1]]), 
                               as.cl_partition(clus_PhenoGraph[[1]]), 
                               as.cl_partition(clus_Rclusterpp[[1]]))

partition_Levine_13dim <- list(as.cl_partition(clus_ClusterX[[2]]), 
                               as.cl_partition(clus_FLOCK[[2]]), 
                               as.cl_partition(clus_flowMeans[[2]]), 
                               as.cl_partition(clus_flowPeaks[[2]]), 
                               as.cl_partition(clus_FlowSOM[[2]]), 
                               as.cl_partition(clus_immunoClust[[2]]), 
                               as.cl_partition(clus_kmeans[[2]]), 
                               as.cl_partition(clus_PhenoGraph[[2]]), 
                               as.cl_partition(clus_Rclusterpp[[2]]), 
                               as.cl_partition(clus_SPADE[[2]]))

partition_Samusik_01   <- list(as.cl_partition(clus_ClusterX[[3]]), 
                               as.cl_partition(clus_DensVM[[3]]), 
                               as.cl_partition(clus_FLOCK[[3]]), 
                               as.cl_partition(clus_flowMeans[[3]]), 
                               as.cl_partition(clus_flowPeaks[[3]]), 
                               as.cl_partition(clus_FlowSOM[[3]]), 
                               as.cl_partition(clus_immunoClust[[3]]), 
                               as.cl_partition(clus_kmeans[[3]]), 
                               as.cl_partition(clus_PhenoGraph[[3]]), 
                               as.cl_partition(clus_Rclusterpp[[3]]), 
                               as.cl_partition(clus_SPADE[[3]]), 
                               as.cl_partition(clus_Xshift[[3]]))

partition_Samusik_all  <- list(as.cl_partition(clus_FLOCK[[4]]), 
                               as.cl_partition(clus_flowPeaks[[4]]), 
                               as.cl_partition(clus_FlowSOM[[4]]), 
                               as.cl_partition(clus_kmeans[[4]]), 
                               as.cl_partition(clus_PhenoGraph[[4]]), 
                               as.cl_partition(clus_SPADE[[4]]))

partition_Nilsson_rare <- list(as.cl_partition(clus_ClusterX[[5]]), 
                               as.cl_partition(clus_DensVM[[5]]), 
                               as.cl_partition(clus_FLOCK[[5]]), 
                               as.cl_partition(clus_flowMeans[[5]]), 
                               as.cl_partition(clus_flowPeaks[[5]]), 
                               as.cl_partition(clus_FlowSOM[[5]]), 
                               as.cl_partition(clus_immunoClust[[5]]), 
                               as.cl_partition(clus_kmeans[[5]]), 
                               as.cl_partition(clus_PhenoGraph[[5]]), 
                               as.cl_partition(clus_Rclusterpp[[5]]), 
                               as.cl_partition(clus_SPADE[[5]]), 
                               as.cl_partition(clus_SWIFT[[5]]), 
                               as.cl_partition(clus_Xshift[[5]]))

partition_Mosmann_rare <- list(as.cl_partition(clus_FLOCK[[6]]), 
                               as.cl_partition(clus_flowMeans[[6]]), 
                               as.cl_partition(clus_flowPeaks[[6]]), 
                               as.cl_partition(clus_FlowSOM[[6]]), 
                               as.cl_partition(clus_immunoClust[[6]]), 
                               as.cl_partition(clus_kmeans[[6]]), 
                               as.cl_partition(clus_PhenoGraph[[6]]), 
                               as.cl_partition(clus_SPADE[[6]]), 
                               as.cl_partition(clus_SWIFT[[6]]))


# create cluster ensembles

partitions <- list(partition_Levine_32dim, 
                   partition_Levine_13dim, 
                   partition_Samusik_01, 
                   partition_Samusik_all, 
                   partition_Nilsson_rare, 
                   partition_Mosmann_rare)

ensembles <- lapply(partitions, function(p) cl_ensemble(list = p))


# calculate consensus clustering (runtime: 35 mins)

consensus <- vector("list", length(ensembles))

for (i in 1:length(consensus)) {
  set.seed(123)
  consensus[[i]] <- cl_consensus(ensembles[[i]])
}


# get class IDs

clus_consensus <- lapply(consensus, cl_class_ids)




####################
### SAVE RESULTS ###
####################

# save cluster labels

datasets <- c("Levine_32dim", "Levine_13dim", "Samusik_01", "Samusik_all", "Nilsson_rare", "Mosmann_rare")
files_out <- paste0("../../results/ensemble/ensemble_labels_", datasets, ".txt")

for (i in 1:length(files_out)) {
  res <- data.frame(label = as.numeric(clus_consensus[[i]]))
  write.table(res, file = files_out[i], row.names = FALSE, quote = FALSE, sep = "\t")
}


# save session information

sink(file = "../../results/ensemble/session_info_ensemble.txt")
sessionInfo()
sink()



