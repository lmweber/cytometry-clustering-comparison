#########################################################################################
# R script to load and calculate results for PhenoGraph
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")



###############################
### load PhenoGraph results ###
###############################

# load PhenoGraph output files

file_PhenoGraph_Levine <- "../results/PhenoGraph/PhenoGraph_results_Levine_2015_marrow_32.fcs"
file_PhenoGraph_Mosmann <- "../results/PhenoGraph/PhenoGraph_results_Mosmann_2014_rare.fcs"

data_PhenoGraph_Levine <- flowCore::exprs(flowCore::read.FCS(file_PhenoGraph_Levine, transformation = FALSE))
data_PhenoGraph_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_PhenoGraph_Mosmann, transformation = FALSE))

head(data_PhenoGraph_Levine)
head(data_PhenoGraph_Mosmann)

dim(data_PhenoGraph_Levine)
dim(data_PhenoGraph_Mosmann)


# extract cluster labels

clus_PhenoGraph_Levine <- data_PhenoGraph_Levine[, grep("PhenoGraph", colnames(data_PhenoGraph_Levine))]
clus_PhenoGraph_Mosmann <- data_PhenoGraph_Mosmann[, grep("PhenoGraph", colnames(data_PhenoGraph_Mosmann))]

length(clus_PhenoGraph_Levine)
length(clus_PhenoGraph_Mosmann)


# contingency tables

table(clus_PhenoGraph_Levine, clus_truth_Levine)
table(clus_PhenoGraph_Mosmann, clus_truth_Mosmann)


# cluster sizes and number of clusters

tbl_PhenoGraph_Levine <- table(clus_PhenoGraph_Levine)
tbl_PhenoGraph_Mosmann <- table(clus_PhenoGraph_Mosmann)

tbl_PhenoGraph_Levine
tbl_PhenoGraph_Mosmann

length(tbl_PhenoGraph_Levine)
length(tbl_PhenoGraph_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_PhenoGraph_Levine <- helper_match_clusters_and_evaluate(clus_PhenoGraph_Levine, clus_truth_Levine)
res_PhenoGraph_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_PhenoGraph_Mosmann, clus_truth_Mosmann)

res_PhenoGraph_Levine
res_PhenoGraph_Mosmann


