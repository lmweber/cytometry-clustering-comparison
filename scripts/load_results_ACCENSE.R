#########################################################################################
# R script to load and calculate results for ACCENSE
#
# Lukas M. Weber, November 2015
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

file_ACCENSE_Levine <- "../results/ACCENSE/accense_output_Levine_2015_marrow_32.csv"
file_ACCENSE_Mosmann <- "../results/ACCENSE/accense_output_Mosmann_2014_rare.csv"

data_ACCENSE_Levine <- read.csv(file_ACCENSE_Levine)
data_ACCENSE_Mosmann <- read.csv(file_ACCENSE_Mosmann)

head(data_ACCENSE_Levine)
head(data_ACCENSE_Mosmann)

dim(data_ACCENSE_Levine)
dim(data_ACCENSE_Mosmann)


# extract truth labels from subsampled data set

clus_truth_Levine_subsampled <- data_ACCENSE_Levine[, "label"]
clus_truth_Mosmann_subsampled <- data_ACCENSE_Mosmann[, "label"]

table(clus_truth_Levine_subsampled)
table(clus_truth_Mosmann_subsampled)  # note: too few cells; subsampling is not suitable for rare cell types

length(clus_truth_Levine_subsampled)
length(clus_truth_Mosmann_subsampled)


# extract ACCENSE cluster labels

clus_ACCENSE_Levine <- data_ACCENSE_Levine[, "population"]
clus_ACCENSE_Mosmann <- data_ACCENSE_Mosmann[, "population"]

length(clus_ACCENSE_Levine)
length(clus_ACCENSE_Mosmann)


# contingency tables

table(clus_ACCENSE_Levine, clus_truth_Levine_subsampled)
table(clus_ACCENSE_Mosmann, clus_truth_Mosmann_subsampled)


# cluster sizes and number of clusters

tbl_ACCENSE_Levine <- table(clus_ACCENSE_Levine)
tbl_ACCENSE_Mosmann <- table(clus_ACCENSE_Mosmann)

tbl_ACCENSE_Levine
tbl_ACCENSE_Mosmann

length(tbl_ACCENSE_Levine)
length(tbl_ACCENSE_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_ACCENSE_Levine <- helper_match_clusters_and_evaluate(clus_ACCENSE_Levine, clus_truth_Levine_subsampled)
res_ACCENSE_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_ACCENSE_Mosmann, clus_truth_Mosmann_subsampled)

res_ACCENSE_Levine
res_ACCENSE_Mosmann


