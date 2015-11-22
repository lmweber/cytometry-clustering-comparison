#########################################################################################
# R script to load and calculate results for SWIFT
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")



##########################
### load SWIFT results ###
##########################

# load SWIFT output files

file_SWIFT_Levine <- "../results/SWIFT/Levine_2015_marrow_32/Levine_2015_marrow_32_notransf.fcs.Cluster_Output.txt"
file_SWIFT_Mosmann <- "../results/SWIFT/Mosmann_2014_rare/Mosmann_2014_rare_notransf.fcs.Cluster_Output.txt"

data_SWIFT_Levine <- read.table(file_SWIFT_Levine, header = TRUE, sep = "\t", comment.char = "")
data_SWIFT_Mosmann <- read.table(file_SWIFT_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_SWIFT_Levine)
head(data_SWIFT_Mosmann)

dim(data_SWIFT_Levine)
dim(data_SWIFT_Mosmann)


# extract cluster labels

clus_SWIFT_Levine <- data_SWIFT_Levine[, "MergeCluster."]
clus_SWIFT_Mosmann <- data_SWIFT_Mosmann[, "MergeCluster."]

length(clus_SWIFT_Levine)
length(clus_SWIFT_Mosmann)


# contingency tables

table(clus_SWIFT_Levine, clus_truth_Levine)
table(clus_SWIFT_Mosmann, clus_truth_Mosmann)


# cluster sizes and number of clusters

tbl_SWIFT_Levine <- table(clus_SWIFT_Levine)
tbl_SWIFT_Mosmann <- table(clus_SWIFT_Mosmann)

tbl_SWIFT_Levine
tbl_SWIFT_Mosmann

length(tbl_SWIFT_Levine)
length(tbl_SWIFT_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_SWIFT_Levine <- helper_match_clusters_and_evaluate(clus_SWIFT_Levine, clus_truth_Levine)
res_SWIFT_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_SWIFT_Mosmann, clus_truth_Mosmann)

res_SWIFT_Levine
res_SWIFT_Mosmann


