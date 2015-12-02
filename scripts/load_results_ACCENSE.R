#########################################################################################
# R script to load and calculate results for ACCENSE
#
# Lukas M. Weber, December 2015
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")



############################
### load ACCENSE results ###
############################

# note ACCENSE uses subsampled data, so truth labels also need to be re-calculated

# load "accense_output.csv" files

file_ACCENSE_Levine_32 <- "../results/ACCENSE/accense_output_Levine_2015_marrow_32.csv"
file_ACCENSE_Levine_13 <- "../results/ACCENSE/accense_output_Levine_2015_marrow_13.csv"
file_ACCENSE_Mosmann <- "../results/ACCENSE/accense_output_Mosmann_2014_rare.csv"

data_ACCENSE_Levine_32 <- read.csv(file_ACCENSE_Levine_32)
data_ACCENSE_Levine_13 <- read.csv(file_ACCENSE_Levine_13)
data_ACCENSE_Mosmann <- read.csv(file_ACCENSE_Mosmann)

head(data_ACCENSE_Levine_32)
head(data_ACCENSE_Levine_13)
head(data_ACCENSE_Mosmann)

dim(data_ACCENSE_Levine_32)
dim(data_ACCENSE_Levine_13)
dim(data_ACCENSE_Mosmann)


# extract truth labels from subsampled data set

clus_truth_Levine_32_subsampled <- data_ACCENSE_Levine_32[, "label"]
clus_truth_Levine_13_subsampled <- data_ACCENSE_Levine_13[, "label"]
clus_truth_Mosmann_subsampled <- data_ACCENSE_Mosmann[, "label"]

table(clus_truth_Levine_32_subsampled)
table(clus_truth_Levine_13_subsampled)
table(clus_truth_Mosmann_subsampled)  # note: too few cells; subsampling is not suitable for rare cell types

length(clus_truth_Levine_32_subsampled)
length(clus_truth_Levine_13_subsampled)
length(clus_truth_Mosmann_subsampled)


# extract ACCENSE cluster labels

clus_ACCENSE_Levine_32 <- data_ACCENSE_Levine_32[, "population"]
clus_ACCENSE_Levine_13 <- data_ACCENSE_Levine_13[, "population"]
clus_ACCENSE_Mosmann <- data_ACCENSE_Mosmann[, "population"]

length(clus_ACCENSE_Levine_32)
length(clus_ACCENSE_Levine_13)
length(clus_ACCENSE_Mosmann)


# contingency tables

table(clus_ACCENSE_Levine_32, clus_truth_Levine_32_subsampled)
table(clus_ACCENSE_Levine_13, clus_truth_Levine_13_subsampled)
table(clus_ACCENSE_Mosmann, clus_truth_Mosmann_subsampled)


# cluster sizes and number of clusters

tbl_ACCENSE_Levine_32 <- table(clus_ACCENSE_Levine_32)
tbl_ACCENSE_Levine_13 <- table(clus_ACCENSE_Levine_13)
tbl_ACCENSE_Mosmann <- table(clus_ACCENSE_Mosmann)

tbl_ACCENSE_Levine_32
tbl_ACCENSE_Levine_13
tbl_ACCENSE_Mosmann

length(tbl_ACCENSE_Levine_32)
length(tbl_ACCENSE_Levine_13)
length(tbl_ACCENSE_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_ACCENSE_Levine_32 <- helper_match_clusters_and_evaluate(clus_ACCENSE_Levine_32, clus_truth_Levine_32_subsampled)
res_ACCENSE_Levine_13 <- helper_match_clusters_and_evaluate(clus_ACCENSE_Levine_13, clus_truth_Levine_13_subsampled)
res_ACCENSE_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_ACCENSE_Mosmann, clus_truth_Mosmann_subsampled)

res_ACCENSE_Levine_32
res_ACCENSE_Levine_13
res_ACCENSE_Mosmann


