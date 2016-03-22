#########################################################################################
# R script to load and calculate results for SWIFT
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")

# results directories
source("load_results_directories.R")



##########################
### load SWIFT results ###
##########################

# load SWIFT output files

file_SWIFT_Levine_32 <- file.path(RES_DIR_SWIFT, "SWIFT/Levine_2015_marrow_32/Levine_2015_marrow_32_notransform.fcs.Cluster_Output.txt")
file_SWIFT_Levine_13 <- file.path(RES_DIR_SWIFT, "SWIFT/Levine_2015_marrow_13/Levine_2015_marrow_13_notransform.fcs.Cluster_Output.txt")
file_SWIFT_Nilsson <- file.path(RES_DIR_SWIFT, "SWIFT/Nilsson_2013_HSC/Nilsson_2013_HSC_notransform.fcs.Cluster_Output.txt")
file_SWIFT_Mosmann <- file.path(RES_DIR_SWIFT, "SWIFT/Mosmann_2014_activ/Mosmann_2014_activ_notransform.fcs.Cluster_Output.txt")

data_SWIFT_Levine_32 <- read.table(file_SWIFT_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_SWIFT_Levine_13 <- read.table(file_SWIFT_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_SWIFT_Nilsson <- read.table(file_SWIFT_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_SWIFT_Mosmann <- read.table(file_SWIFT_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_SWIFT_Levine_32)
head(data_SWIFT_Levine_13)
head(data_SWIFT_Nilsson)
head(data_SWIFT_Mosmann)

dim(data_SWIFT_Levine_32)
dim(data_SWIFT_Levine_13)
dim(data_SWIFT_Nilsson)
dim(data_SWIFT_Mosmann)


# extract cluster labels

clus_SWIFT_Levine_32 <- data_SWIFT_Levine_32[, "MergeCluster."]
clus_SWIFT_Levine_13 <- data_SWIFT_Levine_13[, "MergeCluster."]
clus_SWIFT_Nilsson <- data_SWIFT_Nilsson[, "MergeCluster."]
clus_SWIFT_Mosmann <- data_SWIFT_Mosmann[, "MergeCluster."]

length(clus_SWIFT_Levine_32)
length(clus_SWIFT_Levine_13)
length(clus_SWIFT_Nilsson)
length(clus_SWIFT_Mosmann)


# contingency tables

table(clus_SWIFT_Levine_32, clus_truth_Levine_32)
table(clus_SWIFT_Levine_13, clus_truth_Levine_13)
table(clus_SWIFT_Nilsson, clus_truth_Nilsson)
table(clus_SWIFT_Mosmann, clus_truth_Mosmann)


# cluster sizes and number of clusters

tbl_SWIFT_Levine_32 <- table(clus_SWIFT_Levine_32)
tbl_SWIFT_Levine_13 <- table(clus_SWIFT_Levine_13)
tbl_SWIFT_Nilsson <- table(clus_SWIFT_Nilsson)
tbl_SWIFT_Mosmann <- table(clus_SWIFT_Mosmann)

tbl_SWIFT_Levine_32
tbl_SWIFT_Levine_13
tbl_SWIFT_Nilsson
tbl_SWIFT_Mosmann

length(tbl_SWIFT_Levine_32)
length(tbl_SWIFT_Levine_13)
length(tbl_SWIFT_Nilsson)
length(tbl_SWIFT_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_SWIFT_Levine_32 <- helper_match_clusters_and_evaluate(clus_SWIFT_Levine_32, clus_truth_Levine_32)
res_SWIFT_Levine_13 <- helper_match_clusters_and_evaluate(clus_SWIFT_Levine_13, clus_truth_Levine_13)
res_SWIFT_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_SWIFT_Nilsson, clus_truth_Nilsson)
res_SWIFT_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_SWIFT_Mosmann, clus_truth_Mosmann)

res_SWIFT_Levine_32
res_SWIFT_Levine_13
res_SWIFT_Nilsson
res_SWIFT_Mosmann


