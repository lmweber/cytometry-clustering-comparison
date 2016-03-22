#########################################################################################
# R script to load and calculate results for FLOCK
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
### load FLOCK results ###
##########################

# load FLOCK output files

file_FLOCK_Levine_32 <- file.path(RES_DIR_FLOCK, "FLOCK/flock_results_Levine_2015_marrow_32.txt")
file_FLOCK_Levine_13 <- file.path(RES_DIR_FLOCK, "FLOCK/flock_results_Levine_2015_marrow_13.txt")
file_FLOCK_Nilsson <- file.path(RES_DIR_FLOCK, "FLOCK/flock_results_Nilsson_2013_HSC.txt")
file_FLOCK_Mosmann <- file.path(RES_DIR_FLOCK, "FLOCK/flock_results_Mosmann_2014_activ.txt")

data_FLOCK_Levine_32 <- read.table(file_FLOCK_Levine_32, header = TRUE, sep = "\t")
data_FLOCK_Levine_13 <- read.table(file_FLOCK_Levine_13, header = TRUE, sep = "\t")
data_FLOCK_Nilsson <- read.table(file_FLOCK_Nilsson, header = TRUE, sep = "\t")
data_FLOCK_Mosmann <- read.table(file_FLOCK_Mosmann, header = TRUE, sep = "\t")

head(data_FLOCK_Levine_32)
head(data_FLOCK_Levine_13)
head(data_FLOCK_Nilsson)
head(data_FLOCK_Mosmann)

dim(data_FLOCK_Levine_32)
dim(data_FLOCK_Levine_13)
dim(data_FLOCK_Nilsson)
dim(data_FLOCK_Mosmann)


# extract cluster labels

clus_FLOCK_Levine_32 <- data_FLOCK_Levine_32[, "Population"]
clus_FLOCK_Levine_13 <- data_FLOCK_Levine_13[, "Population"]
clus_FLOCK_Nilsson <- data_FLOCK_Nilsson[, "Population"]
clus_FLOCK_Mosmann <- data_FLOCK_Mosmann[, "Population"]

length(clus_FLOCK_Levine_32)
length(clus_FLOCK_Levine_13)
length(clus_FLOCK_Nilsson)
length(clus_FLOCK_Mosmann)


# contingency tables

table(clus_FLOCK_Levine_32, clus_truth_Levine_32)
table(clus_FLOCK_Levine_13, clus_truth_Levine_13)
table(clus_FLOCK_Nilsson, clus_truth_Nilsson)
table(clus_FLOCK_Mosmann, clus_truth_Mosmann)


# cluster sizes and number of clusters

tbl_FLOCK_Levine_32 <- table(clus_FLOCK_Levine_32)
tbl_FLOCK_Levine_13 <- table(clus_FLOCK_Levine_13)
tbl_FLOCK_Nilsson <- table(clus_FLOCK_Nilsson)
tbl_FLOCK_Mosmann <- table(clus_FLOCK_Mosmann)

tbl_FLOCK_Levine_32
tbl_FLOCK_Levine_13
tbl_FLOCK_Nilsson
tbl_FLOCK_Mosmann

length(tbl_FLOCK_Levine_32)
length(tbl_FLOCK_Levine_13)
length(tbl_FLOCK_Nilsson)
length(tbl_FLOCK_Mosmann)


# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FLOCK_Levine_32 <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine_32, clus_truth_Levine_32)
res_FLOCK_Levine_13 <- helper_match_clusters_and_evaluate(clus_FLOCK_Levine_13, clus_truth_Levine_13)
res_FLOCK_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Nilsson, clus_truth_Nilsson)
res_FLOCK_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FLOCK_Mosmann, clus_truth_Mosmann)

res_FLOCK_Levine_32
res_FLOCK_Levine_13
res_FLOCK_Nilsson
res_FLOCK_Mosmann


