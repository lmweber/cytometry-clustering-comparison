#########################################################################################
# R script to load and calculate results for all other methods
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")



##############################
### load flowMeans results ###
##############################

# load cluster labels

file_flowMeans_Levine <- "../results/flowMeans/flowMeans_labels_Levine_2015_marrow_32.txt"
file_flowMeans_Mosmann <- "../results/flowMeans/flowMeans_labels_Mosmann_2014_rare.txt"

data_flowMeans_Levine <- read.table(file_flowMeans_Levine, header = TRUE, sep = "\t", comment.char = "")
data_flowMeans_Mosmann <- read.table(file_flowMeans_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_flowMeans_Levine)
head(data_flowMeans_Mosmann)

clus_flowMeans_Levine <- data_flowMeans_Levine[, "label"]
clus_flowMeans_Mosmann <- data_flowMeans_Mosmann[, "label"]

length(clus_flowMeans_Levine)
length(clus_flowMeans_Mosmann)

# contingency tables

table(clus_flowMeans_Levine, clus_truth_Levine)
table(clus_flowMeans_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_flowMeans_Levine <- table(clus_flowMeans_Levine)
tbl_flowMeans_Mosmann <- table(clus_flowMeans_Mosmann)

tbl_flowMeans_Levine
tbl_flowMeans_Mosmann

length(tbl_flowMeans_Levine)
length(tbl_flowMeans_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_flowMeans_Levine <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine, clus_truth_Levine)
res_flowMeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Mosmann, clus_truth_Mosmann)

res_flowMeans_Levine
res_flowMeans_Mosmann



############################
### load FlowSOM results ###
############################

# load cluster labels

file_FlowSOM_Levine <- "../results/FlowSOM/FlowSOM_labels_Levine_2015_marrow_32.txt"
file_FlowSOM_Mosmann <- "../results/FlowSOM/FlowSOM_labels_Mosmann_2014_rare.txt"

data_FlowSOM_Levine <- read.table(file_FlowSOM_Levine, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_Mosmann <- read.table(file_FlowSOM_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_FlowSOM_Levine)
head(data_FlowSOM_Mosmann)

clus_FlowSOM_Levine <- data_FlowSOM_Levine[, "label"]
clus_FlowSOM_Mosmann <- data_FlowSOM_Mosmann[, "label"]

length(clus_FlowSOM_Levine)
length(clus_FlowSOM_Mosmann)

# contingency tables

table(clus_FlowSOM_Levine, clus_truth_Levine)
table(clus_FlowSOM_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_FlowSOM_Levine <- table(clus_FlowSOM_Levine)
tbl_FlowSOM_Mosmann <- table(clus_FlowSOM_Mosmann)

tbl_FlowSOM_Levine
tbl_FlowSOM_Mosmann

length(tbl_FlowSOM_Levine)
length(tbl_FlowSOM_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FlowSOM_Levine <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine, clus_truth_Levine)
res_FlowSOM_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Mosmann, clus_truth_Mosmann)

res_FlowSOM_Levine
res_FlowSOM_Mosmann



#################################
### load FlowSOM_meta results ###
#################################

# load cluster labels

file_FlowSOM_meta_Levine <- "../results/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_32.txt"
file_FlowSOM_meta_Mosmann <- "../results/FlowSOM_meta/FlowSOM_meta_labels_Mosmann_2014_rare.txt"

data_FlowSOM_meta_Levine <- read.table(file_FlowSOM_meta_Levine, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_meta_Mosmann <- read.table(file_FlowSOM_meta_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_FlowSOM_meta_Levine)
head(data_FlowSOM_meta_Mosmann)

clus_FlowSOM_meta_Levine <- data_FlowSOM_meta_Levine[, "label"]
clus_FlowSOM_meta_Mosmann <- data_FlowSOM_meta_Mosmann[, "label"]

length(clus_FlowSOM_meta_Levine)
length(clus_FlowSOM_meta_Mosmann)

# contingency tables

table(clus_FlowSOM_meta_Levine, clus_truth_Levine)
table(clus_FlowSOM_meta_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_FlowSOM_meta_Levine <- table(clus_FlowSOM_meta_Levine)
tbl_FlowSOM_meta_Mosmann <- table(clus_FlowSOM_meta_Mosmann)

tbl_FlowSOM_meta_Levine
tbl_FlowSOM_meta_Mosmann

length(tbl_FlowSOM_meta_Levine)
length(tbl_FlowSOM_meta_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FlowSOM_meta_Levine <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine, clus_truth_Levine)
res_FlowSOM_meta_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Mosmann, clus_truth_Mosmann)

res_FlowSOM_meta_Levine
res_FlowSOM_meta_Mosmann



################################
### load immunoClust results ###
################################

# load cluster labels

file_immunoClust_Levine <- "../results/immunoClust/immunoClust_labels_Levine_2015_marrow_32.txt"
file_immunoClust_Mosmann <- "../results/immunoClust/immunoClust_labels_Mosmann_2014_rare.txt"

data_immunoClust_Levine <- read.table(file_immunoClust_Levine, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_Mosmann <- read.table(file_immunoClust_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_immunoClust_Levine)
head(data_immunoClust_Mosmann)

clus_immunoClust_Levine <- data_immunoClust_Levine[, "label"]
clus_immunoClust_Mosmann <- data_immunoClust_Mosmann[, "label"]

length(clus_immunoClust_Levine)
length(clus_immunoClust_Mosmann)

# contingency tables

table(clus_immunoClust_Levine, clus_truth_Levine)
table(clus_immunoClust_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_immunoClust_Levine <- table(clus_immunoClust_Levine)
tbl_immunoClust_Mosmann <- table(clus_immunoClust_Mosmann)

tbl_immunoClust_Levine
tbl_immunoClust_Mosmann

length(tbl_immunoClust_Levine)
length(tbl_immunoClust_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_immunoClust_Levine <- helper_match_clusters_and_evaluate(clus_immunoClust_Levine, clus_truth_Levine)
res_immunoClust_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_Mosmann, clus_truth_Mosmann)

res_immunoClust_Levine
res_immunoClust_Mosmann



####################################
### load immunoClust_all results ###
####################################

# load cluster labels

file_immunoClust_all_Levine <- "../results/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_32.txt"
file_immunoClust_all_Mosmann <- "../results/immunoClust_all/immunoClust_all_labels_Mosmann_2014_rare.txt"

data_immunoClust_all_Levine <- read.table(file_immunoClust_all_Levine, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_all_Mosmann <- read.table(file_immunoClust_all_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_immunoClust_all_Levine)
head(data_immunoClust_all_Mosmann)

clus_immunoClust_all_Levine <- data_immunoClust_all_Levine[, "label"]
clus_immunoClust_all_Mosmann <- data_immunoClust_all_Mosmann[, "label"]

length(clus_immunoClust_all_Levine)
length(clus_immunoClust_all_Mosmann)

# contingency tables

table(clus_immunoClust_all_Levine, clus_truth_Levine)
table(clus_immunoClust_all_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_immunoClust_all_Levine <- table(clus_immunoClust_all_Levine)
tbl_immunoClust_all_Mosmann <- table(clus_immunoClust_all_Mosmann)

tbl_immunoClust_all_Levine
tbl_immunoClust_all_Mosmann

length(tbl_immunoClust_all_Levine)
length(tbl_immunoClust_all_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_immunoClust_all_Levine <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine, clus_truth_Levine)
res_immunoClust_all_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Mosmann, clus_truth_Mosmann)

res_immunoClust_all_Levine
res_immunoClust_all_Mosmann



############################
### load k-means results ###
############################

# load cluster labels

file_kmeans_Levine <- "../results/kmeans/kmeans_labels_Levine_2015_marrow_32.txt"
file_kmeans_Mosmann <- "../results/kmeans/kmeans_labels_Mosmann_2014_rare.txt"

data_kmeans_Levine <- read.table(file_kmeans_Levine, header = TRUE, sep = "\t", comment.char = "")
data_kmeans_Mosmann <- read.table(file_kmeans_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_kmeans_Levine)
head(data_kmeans_Mosmann)

clus_kmeans_Levine <- data_kmeans_Levine[, "label"]
clus_kmeans_Mosmann <- data_kmeans_Mosmann[, "label"]

length(clus_kmeans_Levine)
length(clus_kmeans_Mosmann)

# contingency tables

table(clus_kmeans_Levine, clus_truth_Levine)
table(clus_kmeans_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_kmeans_Levine <- table(clus_kmeans_Levine)
tbl_kmeans_Mosmann <- table(clus_kmeans_Mosmann)

tbl_kmeans_Levine
tbl_kmeans_Mosmann

length(tbl_kmeans_Levine)
length(tbl_kmeans_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_kmeans_Levine <- helper_match_clusters_and_evaluate(clus_kmeans_Levine, clus_truth_Levine)
res_kmeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_kmeans_Mosmann, clus_truth_Mosmann)

res_kmeans_Levine
res_kmeans_Mosmann



###############################
### load Rclusterpp results ###
###############################

# load cluster labels

file_Rclusterpp_Levine <- "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_32.txt"
# file_Rclusterpp_Mosmann <- "../results/Rclusterpp/Rclusterpp_labels_Mosmann_2014_rare.txt"

data_Rclusterpp_Levine <- read.table(file_Rclusterpp_Levine, header = TRUE, sep = "\t", comment.char = "")
# data_Rclusterpp_Mosmann <- read.table(file_Rclusterpp_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_Rclusterpp_Levine)
# head(data_Rclusterpp_Mosmann)

clus_Rclusterpp_Levine <- data_Rclusterpp_Levine[, "label"]
# clus_Rclusterpp_Mosmann <- data_Rclusterpp_Mosmann[, "label"]

length(clus_Rclusterpp_Levine)
# length(clus_Rclusterpp_Mosmann)

# contingency tables

table(clus_Rclusterpp_Levine, clus_truth_Levine)
# table(clus_Rclusterpp_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_Rclusterpp_Levine <- table(clus_Rclusterpp_Levine)
# tbl_Rclusterpp_Mosmann <- table(clus_Rclusterpp_Mosmann)

tbl_Rclusterpp_Levine
# tbl_Rclusterpp_Mosmann

length(tbl_Rclusterpp_Levine)
# length(tbl_Rclusterpp_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_Rclusterpp_Levine <- helper_match_clusters_and_evaluate(clus_Rclusterpp_Levine, clus_truth_Levine)
# res_Rclusterpp_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_Rclusterpp_Mosmann, clus_truth_Mosmann)

res_Rclusterpp_Levine
# res_Rclusterpp_Mosmann



################################
### load SamSPECTRAL results ###
################################

# load cluster labels

file_SamSPECTRAL_Levine <- "../results/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_32.txt"
file_SamSPECTRAL_Mosmann <- "../results/SamSPECTRAL/SamSPECTRAL_labels_Mosmann_2014_rare.txt"

data_SamSPECTRAL_Levine <- read.table(file_SamSPECTRAL_Levine, header = TRUE, sep = "\t", comment.char = "")
data_SamSPECTRAL_Mosmann <- read.table(file_SamSPECTRAL_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_SamSPECTRAL_Levine)
head(data_SamSPECTRAL_Mosmann)

clus_SamSPECTRAL_Levine <- data_SamSPECTRAL_Levine[, "label"]
clus_SamSPECTRAL_Mosmann <- data_SamSPECTRAL_Mosmann[, "label"]

length(clus_SamSPECTRAL_Levine)
length(clus_SamSPECTRAL_Mosmann)

# contingency tables

table(clus_SamSPECTRAL_Levine, clus_truth_Levine)
table(clus_SamSPECTRAL_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_SamSPECTRAL_Levine <- table(clus_SamSPECTRAL_Levine)
tbl_SamSPECTRAL_Mosmann <- table(clus_SamSPECTRAL_Mosmann)

tbl_SamSPECTRAL_Levine
tbl_SamSPECTRAL_Mosmann

length(tbl_SamSPECTRAL_Levine)
length(tbl_SamSPECTRAL_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_SamSPECTRAL_Levine <- helper_match_clusters_and_evaluate(clus_SamSPECTRAL_Levine, clus_truth_Levine)
res_SamSPECTRAL_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_SamSPECTRAL_Mosmann, clus_truth_Mosmann)

res_SamSPECTRAL_Levine
res_SamSPECTRAL_Mosmann


