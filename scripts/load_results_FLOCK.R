#########################################################################################
# R script to load and calculate results for FLOCK
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
### load FLOCK results ###
##########################

# load FLOCK output files

file_FLOCK_Levine <- "../results/FLOCK/flock_results_Levine_2015_marrow_32.txt"
file_FLOCK_Mosmann <- "../results/FLOCK/flock_results_Mosmann_2014_rare.txt"

data_FLOCK_Levine <- read.table(file_FLOCK_Levine, header = TRUE, sep = "\t")
data_FLOCK_Mosmann <- read.table(file_FLOCK_Mosmann, header = TRUE, sep = "\t")

head(data_FLOCK_Levine)
head(data_FLOCK_Mosmann)

dim(data_FLOCK_Levine)
dim(data_FLOCK_Mosmann)


# extract cluster labels

clus_FLOCK_Levine <- data_FLOCK_Levine[, "Population"]
clus_FLOCK_Mosmann <- data_FLOCK_Mosmann[, "Population"]

length(clus_FLOCK_Levine)
length(clus_FLOCK_Mosmann)


# contingency tables

table(clus_FLOCK_Levine, clus_truth_Levine)
table(clus_FLOCK_Mosmann, clus_truth_Mosmann)


# cluster sizes and number of clusters

tbl_FLOCK_Levine <- table(clus_FLOCK_Levine)
tbl_FLOCK_Mosmann <- table(clus_FLOCK_Mosmann)

tbl_FLOCK_Levine
tbl_FLOCK_Mosmann

length(tbl_FLOCK_Levine)
length(tbl_FLOCK_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FLOCK_Levine <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine, clus_truth_Levine)
res_FLOCK_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Mosmann, clus_truth_Mosmann)

res_FLOCK_Levine
res_FLOCK_Mosmann


