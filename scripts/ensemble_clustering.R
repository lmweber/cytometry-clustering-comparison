#########################################################################################
# R script to run ensemble clustering
#
# Lukas M. Weber, February 2016
#########################################################################################


library(clue)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# load clustering results from previous steps
source("load_results_ACCENSE.R")
source("load_results_DensVM.R")
source("load_results_FLOCK.R")
source("load_results_PhenoGraph.R")
source("load_results_SWIFT.R")
source("load_results_truth.R")
source("load_results_all_other_methods.R")




###########################
### ENSEMBLE CLUSTERING ###
###########################

# using top 4 methods by mean F1 score for each data set (excluding SamSPECTRAL since it
# gives an error)

# Note: Can't use methods that require subsampled data (ACCENSE, DensVM), since require 
# same data for all methods. For data sets with multiple clusters (Levine_32, Levine_13),
# it is also difficult to use methods that return too many clusters (FlowSOM, 
# immunoClust, immunoClust_all, SWIFT).


# create partitions in format required by CLUE

partition_Levine_32 <- list(as.cl_partition(clus_FlowSOM_meta_Levine_32), 
                            as.cl_partition(clus_flowMeans_Levine_32), 
                            as.cl_partition(clus_FLOCK_Levine_32), 
                            as.cl_partition(clus_PhenoGraph_Levine_32))

partition_Levine_13 <- list(as.cl_partition(clus_FlowSOM_meta_Levine_13), 
                            as.cl_partition(clus_PhenoGraph_Levine_13), 
                            as.cl_partition(clus_flowMeans_Levine_13), 
                            as.cl_partition(clus_FLOCK_Levine_13))

partition_Nilsson <- list(as.cl_partition(clus_flowMeans_Nilsson), 
                          as.cl_partition(clus_FlowSOM_meta_Nilsson), 
                          as.cl_partition(clus_FlowSOM_Nilsson), 
                          as.cl_partition(clus_kmeans_Nilsson))

# excluding SamSPECTRAL since it gives an error
partition_Mosmann <- list(as.cl_partition(clus_SWIFT_Mosmann), 
                          #as.cl_partition(clus_SamSPECTRAL_Mosmann),  ## gives error
                          as.cl_partition(clus_PhenoGraph_Mosmann), 
                          as.cl_partition(clus_FlowSOM_meta_Mosmann))


# create cluster ensembles

ens_Levine_32 <- cl_ensemble(list = partition_Levine_32)
ens_Levine_13 <- cl_ensemble(list = partition_Levine_13)
ens_Nilsson <- cl_ensemble(list = partition_Nilsson)
ens_Mosmann <- cl_ensemble(list = partition_Mosmann)


# calculate consensus clustering

set.seed(1234)
system.time(
  cons_Levine_32 <- cl_consensus(ens_Levine_32)  # runtime: ~5 sec
)

set.seed(1234)
system.time(
  cons_Levine_13 <- cl_consensus(ens_Levine_13)  # runtime: ~6 sec
)

set.seed(1234)
system.time(
  cons_Nilsson <- cl_consensus(ens_Nilsson)  # runtime: ~19 sec
)

set.seed(1234)
system.time(
  cons_Mosmann <- cl_consensus(ens_Mosmann)  # runtime: ~4 min
)


# get class IDs

cons_clus_Levine_32 <- cl_class_ids(cons_Levine_32)
cons_clus_Levine_13 <- cl_class_ids(cons_Levine_13)
cons_clus_Nilsson <- cl_class_ids(cons_Nilsson)
cons_clus_Mosmann <- cl_class_ids(cons_Mosmann)


# cluster frequency tables

table(cons_clus_Levine_32)
table(cons_clus_Levine_13)
table(cons_clus_Nilsson)
table(cons_clus_Mosmann)

length(table(cons_clus_Nilsson))  # 94 clusters
length(table(cons_clus_Mosmann))  # 147 clusters


# confusion matrices

table(ensemble = cons_clus_Levine_32, true = clus_truth_Levine_32)
table(ensemble = cons_clus_Levine_13, true = clus_truth_Levine_13)
table(ensemble = cons_clus_Nilsson, true = clus_truth_Nilsson)
table(ensemble = cons_clus_Mosmann, true = clus_truth_Mosmann)




###############
### RESULTS ###
###############

# for data sets with multiple populations of interest (Levine_32, Levine_13)

# calculate F1 score, precision, recall

res_ensemble_Levine_32 <- helper_match_clusters_and_evaluate(cons_clus_Levine_32, clus_truth_Levine_32)
res_ensemble_Levine_13 <- helper_match_clusters_and_evaluate(cons_clus_Levine_13, clus_truth_Levine_13)

res_ensemble_Levine_32
res_ensemble_Levine_13

mean_F1_ensemble_Levine_32 <- mean(res_ensemble_Levine_32$F1)
mean_F1_ensemble_Levine_13 <- mean(res_ensemble_Levine_13$F1)

mean_F1_ensemble_Levine_32
mean_F1_ensemble_Levine_13


# for data sets with a single rare cell population of interest (Nilsson, Mosmann)

# calculate F1 score, precision, recall

res_ensemble_Nilsson <- helper_match_one_rare_cluster_and_evaluate(cons_clus_Nilsson, clus_truth_Nilsson)
res_ensemble_Mosmann <- helper_match_one_rare_cluster_and_evaluate(cons_clus_Mosmann, clus_truth_Mosmann)

as.data.frame(res_ensemble_Nilsson)
as.data.frame(res_ensemble_Mosmann)




###################################
### SESSION INFO AND RDATA FILE ###
###################################

# save session information

sink(file = "../results/session_info/ensemble_clustering_session_info.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/ensemble_clustering_results.RData")



