#########################################################################################
# R script to load and calculate results for DensVM
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")



###########################
### load DensVM results ###
###########################

# note DensVM uses subsampled data, so truth labels also need to be re-calculated

# load DensVM output files

file_DensVM_Levine <- "../results/DensVM/Levine_2015_marrow_32/cytofkit_analysis_tsne_cluster.txt"
file_DensVM_Mosmann <- "../results/DensVM/Mosmann_2014_rare/cytofkit_analysis_tsne_cluster.txt"

data_DensVM_Levine <- read.table(file_DensVM_Levine, header = TRUE, sep = "\t")
data_DensVM_Mosmann <- read.table(file_DensVM_Mosmann, header = TRUE, sep = "\t")

head(data_DensVM_Levine)
head(data_DensVM_Mosmann)

dim(data_DensVM_Levine)
dim(data_DensVM_Mosmann)


# extract truth labels from subsampled data set

file_truth_Levine_subsampled <- "../results/DensVM/Levine_2015_marrow_32/cytofkit_analysis_analyzedFCS/Levine_2015_marrow_32_notransf.fcs"
file_truth_Mosmann_subsampled <- "../results/DensVM/Mosmann_2014_rare/cytofkit_analysis_analyzedFCS/Mosmann_2014_rare_notransf.fcs"

data_truth_Levine_subsampled <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_subsampled, transformation = FALSE))
data_truth_Mosmann_subsampled <- flowCore::exprs(flowCore::read.FCS(file_truth_Mosmann_subsampled, transformation = FALSE))

dim(data_truth_Levine_subsampled)
dim(data_truth_Mosmann_subsampled)

clus_truth_Levine_subsampled <- data_truth_Levine_subsampled[, "label"]
clus_truth_Mosmann_subsampled <- data_truth_Mosmann_subsampled[, "label"]

table(clus_truth_Levine_subsampled)
table(clus_truth_Mosmann_subsampled)  # note: too few cells; subsampling is not suitable for rare cell types

length(clus_truth_Levine_subsampled)
length(clus_truth_Mosmann_subsampled)


# extract DensVM cluster labels

clus_DensVM_Levine <- data_DensVM_Levine[, "cluster"]
clus_DensVM_Mosmann <- data_DensVM_Mosmann[, "cluster"]

length(clus_DensVM_Levine)
length(clus_DensVM_Mosmann)


# contingency tables

table(clus_DensVM_Levine, clus_truth_Levine_subsampled)
table(clus_DensVM_Mosmann, clus_truth_Mosmann_subsampled)


# cluster sizes and number of clusters

tbl_DensVM_Levine <- table(clus_DensVM_Levine)
tbl_DensVM_Mosmann <- table(clus_DensVM_Mosmann)

tbl_DensVM_Levine
tbl_DensVM_Mosmann

length(tbl_DensVM_Levine)
length(tbl_DensVM_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_DensVM_Levine <- helper_match_clusters_and_evaluate(clus_DensVM_Levine, clus_truth_Levine_subsampled)
res_DensVM_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_DensVM_Mosmann, clus_truth_Mosmann_subsampled)

res_DensVM_Levine
res_DensVM_Mosmann


