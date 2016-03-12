#########################################################################################
# R script to load and calculate results for Rclusterpp
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")



###############################
### load Rclusterpp results ###
###############################

# note Rclusterpp uses subsampled data, so truth labels also need to be re-calculated

# load true population labels from subsampled data

file_truth_Levine_32_sub <- "../results/Rclusterpp/Levine_2015_marrow_32_sub.fcs"
file_truth_Levine_13_sub <- "../results/Rclusterpp/Levine_2015_marrow_13_sub.fcs"
file_truth_Nilsson_sub <- "../results/Rclusterpp/Nilsson_2013_HSC_sub.fcs"
file_truth_Mosmann_sub <- "../results/Rclusterpp/Mosmann_2014_activ_sub.fcs"

data_truth_Levine_32_sub <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_32_sub, transformation = FALSE))
data_truth_Levine_13_sub <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_13_sub, transformation = FALSE))
data_truth_Nilsson_sub <- flowCore::exprs(flowCore::read.FCS(file_truth_Nilsson_sub, transformation = FALSE))
data_truth_Mosmann_sub <- flowCore::exprs(flowCore::read.FCS(file_truth_Mosmann_sub, transformation = FALSE))

dim(data_truth_Levine_32_sub)
dim(data_truth_Levine_13_sub)
dim(data_truth_Nilsson_sub)
dim(data_truth_Mosmann_sub)

clus_truth_Levine_32_sub <- data_truth_Levine_32_sub[, "label"]
clus_truth_Levine_13_sub <- data_truth_Levine_13_sub[, "label"]
clus_truth_Nilsson_sub <- data_truth_Nilsson_sub[, "label"]
clus_truth_Mosmann_sub <- data_truth_Mosmann_sub[, "label"]

table(clus_truth_Levine_32_sub)
table(clus_truth_Levine_13_sub)
table(clus_truth_Nilsson_sub)
table(clus_truth_Mosmann_sub)

length(clus_truth_Levine_32_sub)
length(clus_truth_Levine_13_sub)
length(clus_truth_Nilsson_sub)
length(clus_truth_Mosmann_sub)


# load Rclusterpp cluster labels

file_Rclusterpp_Levine_32 <- "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_32.txt"
file_Rclusterpp_Levine_13 <- "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_13.txt"
file_Rclusterpp_Nilsson <- "../results/Rclusterpp/Rclusterpp_labels_Nilsson_2013_HSC.txt"
file_Rclusterpp_Mosmann <- "../results/Rclusterpp/Rclusterpp_labels_Mosmann_2014_activ.txt"

data_Rclusterpp_Levine_32 <- read.table(file_Rclusterpp_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_Rclusterpp_Levine_13 <- read.table(file_Rclusterpp_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_Rclusterpp_Nilsson <- read.table(file_Rclusterpp_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_Rclusterpp_Mosmann <- read.table(file_Rclusterpp_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_Rclusterpp_Levine_32)
head(data_Rclusterpp_Levine_13)
head(data_Rclusterpp_Nilsson)
head(data_Rclusterpp_Mosmann)

clus_Rclusterpp_Levine_32 <- data_Rclusterpp_Levine_32[, "label"]
clus_Rclusterpp_Levine_13 <- data_Rclusterpp_Levine_13[, "label"]
clus_Rclusterpp_Nilsson <- data_Rclusterpp_Nilsson[, "label"]
clus_Rclusterpp_Mosmann <- data_Rclusterpp_Mosmann[, "label"]

length(clus_Rclusterpp_Levine_32)
length(clus_Rclusterpp_Levine_13)
length(clus_Rclusterpp_Nilsson)
length(clus_Rclusterpp_Mosmann)


# contingency tables

table(clus_Rclusterpp_Levine_32, clus_truth_Levine_32_sub)
table(clus_Rclusterpp_Levine_13, clus_truth_Levine_13_sub)
table(clus_Rclusterpp_Nilsson, clus_truth_Nilsson_sub)
table(clus_Rclusterpp_Mosmann, clus_truth_Mosmann_sub)


# cluster sizes and number of clusters

tbl_Rclusterpp_Levine_32 <- table(clus_Rclusterpp_Levine_32)
tbl_Rclusterpp_Levine_13 <- table(clus_Rclusterpp_Levine_13)
tbl_Rclusterpp_Nilsson <- table(clus_Rclusterpp_Nilsson)
tbl_Rclusterpp_Mosmann <- table(clus_Rclusterpp_Mosmann)

tbl_Rclusterpp_Levine_32
tbl_Rclusterpp_Levine_13
tbl_Rclusterpp_Nilsson
tbl_Rclusterpp_Mosmann

length(tbl_Rclusterpp_Levine_32)
length(tbl_Rclusterpp_Levine_13)
length(tbl_Rclusterpp_Nilsson)
length(tbl_Rclusterpp_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_Rclusterpp_Levine_32 <- helper_match_clusters_and_evaluate(clus_Rclusterpp_Levine_32, clus_truth_Levine_32_sub)
res_Rclusterpp_Levine_13 <- helper_match_clusters_and_evaluate(clus_Rclusterpp_Levine_13, clus_truth_Levine_13_sub)
res_Rclusterpp_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_Rclusterpp_Nilsson, clus_truth_Nilsson_sub)
res_Rclusterpp_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_Rclusterpp_Mosmann, clus_truth_Mosmann_sub)

res_Rclusterpp_Levine_32
res_Rclusterpp_Levine_13
res_Rclusterpp_Nilsson
res_Rclusterpp_Mosmann


