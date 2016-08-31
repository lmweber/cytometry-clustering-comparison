#########################################################################################
# R script to run ensemble clustering
#
# Lukas Weber, August 2016
#########################################################################################


library(clue)

# load clustering results from previous steps
source("load_results_FLOCK.R")
source("load_results_PhenoGraph.R")
source("load_results_SWIFT.R")
source("load_results_truth.R")
source("load_results_all_other_methods.R")




###########################
### ENSEMBLE CLUSTERING ###
###########################

# calculate using the top 5 methods for each data set, excluding the following:
# - methods that use subsampled data (ACCENSE, DensVM, Rclusterpp) and methods that 
# remove outliers (ACCENSE, immunoClust, SamSPECTRAL), since consensus clustering
# requires same data for each method
# - methods that give large numbers of clusters (FlowSOM, immunoClust, immunoClust_all,
# SWIFT), since this greatly increases runtime; except for the Nilsson data set, which is
# small enough that runtime remains fast


# create partitions in format required by CLUE

partition_Levine_32 <- list(as.cl_partition(clus_FlowSOM_meta_Levine_32), 
                            as.cl_partition(clus_flowMeans_Levine_32), 
                            as.cl_partition(clus_FLOCK_Levine_32), 
                            as.cl_partition(clus_PhenoGraph_Levine_32), 
                            as.cl_partition(clus_kmeans_Levine_32))

partition_Levine_13 <- list(as.cl_partition(clus_FlowSOM_meta_Levine_13), 
                            as.cl_partition(clus_PhenoGraph_Levine_13), 
                            as.cl_partition(clus_flowMeans_Levine_13), 
                            as.cl_partition(clus_kmeans_Levine_13), 
                            as.cl_partition(clus_FLOCK_Levine_13))

partition_Nilsson <- list(as.cl_partition(clus_kmeans_Nilsson), 
                          as.cl_partition(clus_FlowSOM_meta_Nilsson), 
                          as.cl_partition(clus_FlowSOM_Nilsson), 
                          as.cl_partition(clus_flowMeans_Nilsson), 
                          as.cl_partition(clus_SWIFT_Nilsson))

partition_Mosmann <- list(as.cl_partition(clus_FlowSOM_meta_Mosmann), 
                          as.cl_partition(clus_PhenoGraph_Mosmann), 
                          as.cl_partition(clus_flowMeans_Mosmann), 
                          as.cl_partition(clus_FLOCK_Mosmann), 
                          as.cl_partition(clus_kmeans_Mosmann))


# create cluster ensembles

ens_Levine_32 <- cl_ensemble(list = partition_Levine_32)
ens_Levine_13 <- cl_ensemble(list = partition_Levine_13)
ens_Nilsson <- cl_ensemble(list = partition_Nilsson)
ens_Mosmann <- cl_ensemble(list = partition_Mosmann)


# calculate consensus clustering

set.seed(1234)
system.time(
  cons_Levine_32 <- cl_consensus(ens_Levine_32)  # runtime: ~1 min
)

set.seed(1234)
system.time(
  cons_Levine_13 <- cl_consensus(ens_Levine_13)  # runtime: ~20 sec
)

set.seed(1234)
system.time(
  cons_Nilsson <- cl_consensus(ens_Nilsson)  # runtime: ~20 sec
)

set.seed(1234)
system.time(
  cons_Mosmann <- cl_consensus(ens_Mosmann)  # runtime: ~1 min
)


# get class IDs

cons_clus_Levine_32 <- cl_class_ids(cons_Levine_32)
cons_clus_Levine_13 <- cl_class_ids(cons_Levine_13)
cons_clus_Nilsson <- cl_class_ids(cons_Nilsson)
cons_clus_Mosmann <- cl_class_ids(cons_Mosmann)




####################
### SAVE RESULTS ###
####################

# save cluster labels

res_ensemble_Levine_32 <- data.frame(label = as.numeric(cons_clus_Levine_32))
res_ensemble_Levine_13 <- data.frame(label = as.numeric(cons_clus_Levine_13))
res_ensemble_Nilsson <- data.frame(label = as.numeric(cons_clus_Nilsson))
res_ensemble_Mosmann <- data.frame(label = as.numeric(cons_clus_Mosmann))

write.table(res_ensemble_Levine_32, 
            file = "../results_ensemble/ensemble_clustering/ensemble_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_ensemble_Levine_13, 
            file = "../results_ensemble/ensemble_clustering/ensemble_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_ensemble_Nilsson, 
            file = "../results_ensemble/ensemble_clustering/ensemble_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_ensemble_Mosmann, 
            file = "../results_ensemble/ensemble_clustering/ensemble_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_ensemble/session_info/session_info_ensemble_clustering.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_ensemble/RData_files/results_ensemble_clustering.RData")



