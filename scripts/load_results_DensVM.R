#########################################################################################
# R script to load and calculate results for DensVM
#
# Lukas M. Weber, December 2015
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

file_DensVM_Levine_32 <- "../results/DensVM/Levine_2015_marrow_32/cytofkit_analysis_tsne_cluster.txt"
file_DensVM_Levine_13 <- "../results/DensVM/Levine_2015_marrow_13/cytofkit_analysis_tsne_cluster.txt"
file_DensVM_Nilsson <- "../results/DensVM/Nilsson_2013_HSC/cytofkit_analysis_tsne_cluster.txt"
file_DensVM_Mosmann <- "../results/DensVM/Mosmann_2014_activ/cytofkit_analysis_tsne_cluster.txt"

data_DensVM_Levine_32 <- read.table(file_DensVM_Levine_32, header = TRUE, sep = "\t")
data_DensVM_Levine_13 <- read.table(file_DensVM_Levine_13, header = TRUE, sep = "\t")
data_DensVM_Nilsson <- read.table(file_DensVM_Nilsson, header = TRUE, sep = "\t")
data_DensVM_Mosmann <- read.table(file_DensVM_Mosmann, header = TRUE, sep = "\t")

head(data_DensVM_Levine_32)
head(data_DensVM_Levine_13)
head(data_DensVM_Nilsson)
head(data_DensVM_Mosmann)

dim(data_DensVM_Levine_32)
dim(data_DensVM_Levine_13)
dim(data_DensVM_Nilsson)
dim(data_DensVM_Mosmann)


# extract truth labels from subsampled data set

file_truth_Levine_32_subsampled <- "../results/DensVM/Levine_2015_marrow_32/cytofkit_analysis_analyzedFCS/Levine_2015_marrow_32_notransform.fcs"
file_truth_Levine_13_subsampled <- "../results/DensVM/Levine_2015_marrow_13/cytofkit_analysis_analyzedFCS/Levine_2015_marrow_13_notransform.fcs"
file_truth_Nilsson_subsampled <- "../results/DensVM/Nilsson_2013_HSC/cytofkit_analysis_analyzedFCS/Nilsson_2013_HSC_notransform.fcs"
file_truth_Mosmann_subsampled <- "../results/DensVM/Mosmann_2014_activ/cytofkit_analysis_analyzedFCS/Mosmann_2014_activ_notransform.fcs"

data_truth_Levine_32_subsampled <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_32_subsampled, transformation = FALSE))
data_truth_Levine_13_subsampled <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_13_subsampled, transformation = FALSE))
data_truth_Nilsson_subsampled <- flowCore::exprs(flowCore::read.FCS(file_truth_Nilsson_subsampled, transformation = FALSE))
data_truth_Mosmann_subsampled <- flowCore::exprs(flowCore::read.FCS(file_truth_Mosmann_subsampled, transformation = FALSE))

dim(data_truth_Levine_32_subsampled)
dim(data_truth_Levine_13_subsampled)
dim(data_truth_Nilsson_subsampled)
dim(data_truth_Mosmann_subsampled)

clus_truth_Levine_32_subsampled <- data_truth_Levine_32_subsampled[, "label"]
clus_truth_Levine_13_subsampled <- data_truth_Levine_13_subsampled[, "label"]
clus_truth_Nilsson_subsampled <- data_truth_Nilsson_subsampled[, "label"]
clus_truth_Mosmann_subsampled <- data_truth_Mosmann_subsampled[, "label"]

table(clus_truth_Levine_32_subsampled)
table(clus_truth_Levine_13_subsampled)
table(clus_truth_Nilsson_subsampled)
table(clus_truth_Mosmann_subsampled)  # note: too few cells; subsampling is not suitable for rare cell types

length(clus_truth_Levine_32_subsampled)
length(clus_truth_Levine_13_subsampled)
length(clus_truth_Nilsson_subsampled)
length(clus_truth_Mosmann_subsampled)


# extract DensVM cluster labels

clus_DensVM_Levine_32 <- data_DensVM_Levine_32[, "cluster"]
clus_DensVM_Levine_13 <- data_DensVM_Levine_13[, "cluster"]
clus_DensVM_Nilsson <- data_DensVM_Nilsson[, "cluster"]
clus_DensVM_Mosmann <- data_DensVM_Mosmann[, "cluster"]

length(clus_DensVM_Levine_32)
length(clus_DensVM_Levine_13)
length(clus_DensVM_Nilsson)
length(clus_DensVM_Mosmann)


# contingency tables

table(clus_DensVM_Levine_32, clus_truth_Levine_32_subsampled)
table(clus_DensVM_Levine_13, clus_truth_Levine_13_subsampled)
table(clus_DensVM_Nilsson, clus_truth_Nilsson_subsampled)
table(clus_DensVM_Mosmann, clus_truth_Mosmann_subsampled)


# cluster sizes and number of clusters

tbl_DensVM_Levine_32 <- table(clus_DensVM_Levine_32)
tbl_DensVM_Levine_13 <- table(clus_DensVM_Levine_13)
tbl_DensVM_Nilsson <- table(clus_DensVM_Nilsson)
tbl_DensVM_Mosmann <- table(clus_DensVM_Mosmann)

tbl_DensVM_Levine_32
tbl_DensVM_Levine_13
tbl_DensVM_Nilsson
tbl_DensVM_Mosmann

length(tbl_DensVM_Levine_32)
length(tbl_DensVM_Levine_13)
length(tbl_DensVM_Nilsson)
length(tbl_DensVM_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_DensVM_Levine_32 <- helper_match_clusters_and_evaluate(clus_DensVM_Levine_32, clus_truth_Levine_32_subsampled)
res_DensVM_Levine_13 <- helper_match_clusters_and_evaluate(clus_DensVM_Levine_13, clus_truth_Levine_13_subsampled)
res_DensVM_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_DensVM_Nilsson, clus_truth_Nilsson_subsampled)
res_DensVM_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_DensVM_Mosmann, clus_truth_Mosmann_subsampled)

res_DensVM_Levine_32
res_DensVM_Levine_13
res_DensVM_Nilsson
res_DensVM_Mosmann


