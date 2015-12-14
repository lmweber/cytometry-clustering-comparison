#########################################################################################
# R script to load and calculate results for all other methods
#
# Lukas M. Weber, December 2015
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

file_flowMeans_Levine_32 <- "../results/flowMeans/flowMeans_labels_Levine_2015_marrow_32.txt"
file_flowMeans_Levine_13 <- "../results/flowMeans/flowMeans_labels_Levine_2015_marrow_13.txt"
file_flowMeans_Nilsson <- "../results/flowMeans/flowMeans_labels_Nilsson_2013_HSC.txt"
file_flowMeans_Mosmann <- "../results/flowMeans/flowMeans_labels_Mosmann_2014_activ.txt"

data_flowMeans_Levine_32 <- read.table(file_flowMeans_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_flowMeans_Levine_13 <- read.table(file_flowMeans_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_flowMeans_Nilsson <- read.table(file_flowMeans_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_flowMeans_Mosmann <- read.table(file_flowMeans_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_flowMeans_Levine_32)
head(data_flowMeans_Levine_13)
head(data_flowMeans_Nilsson)
head(data_flowMeans_Mosmann)

clus_flowMeans_Levine_32 <- data_flowMeans_Levine_32[, "label"]
clus_flowMeans_Levine_13 <- data_flowMeans_Levine_13[, "label"]
clus_flowMeans_Nilsson <- data_flowMeans_Nilsson[, "label"]
clus_flowMeans_Mosmann <- data_flowMeans_Mosmann[, "label"]

length(clus_flowMeans_Levine_32)
length(clus_flowMeans_Levine_13)
length(clus_flowMeans_Nilsson)
length(clus_flowMeans_Mosmann)

# contingency tables

table(clus_flowMeans_Levine_32, clus_truth_Levine_32)
table(clus_flowMeans_Levine_13, clus_truth_Levine_13)
table(clus_flowMeans_Nilsson, clus_truth_Nilsson)
table(clus_flowMeans_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_flowMeans_Levine_32 <- table(clus_flowMeans_Levine_32)
tbl_flowMeans_Levine_13 <- table(clus_flowMeans_Levine_13)
tbl_flowMeans_Nilsson <- table(clus_flowMeans_Nilsson)
tbl_flowMeans_Mosmann <- table(clus_flowMeans_Mosmann)

tbl_flowMeans_Levine_32
tbl_flowMeans_Levine_13
tbl_flowMeans_Nilsson
tbl_flowMeans_Mosmann

length(tbl_flowMeans_Levine_32)
length(tbl_flowMeans_Levine_13)
length(tbl_flowMeans_Nilsson)
length(tbl_flowMeans_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_flowMeans_Levine_32 <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine_32, clus_truth_Levine_32)
res_flowMeans_Levine_13 <- helper_match_clusters_and_evaluate(clus_flowMeans_Levine_13, clus_truth_Levine_13)
res_flowMeans_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Nilsson, clus_truth_Nilsson)
res_flowMeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_flowMeans_Mosmann, clus_truth_Mosmann)

res_flowMeans_Levine_32
res_flowMeans_Levine_13
res_flowMeans_Nilsson
res_flowMeans_Mosmann




############################
### load FlowSOM results ###
############################

# load cluster labels

file_FlowSOM_Levine_32 <- "../results/FlowSOM/FlowSOM_labels_Levine_2015_marrow_32.txt"
file_FlowSOM_Levine_13 <- "../results/FlowSOM/FlowSOM_labels_Levine_2015_marrow_13.txt"
file_FlowSOM_Nilsson <- "../results/FlowSOM/FlowSOM_labels_Nilsson_2013_HSC.txt"
file_FlowSOM_Mosmann <- "../results/FlowSOM/FlowSOM_labels_Mosmann_2014_activ.txt"

data_FlowSOM_Levine_32 <- read.table(file_FlowSOM_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_Levine_13 <- read.table(file_FlowSOM_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_Nilsson <- read.table(file_FlowSOM_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_Mosmann <- read.table(file_FlowSOM_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_FlowSOM_Levine_32)
head(data_FlowSOM_Levine_13)
head(data_FlowSOM_Nilsson)
head(data_FlowSOM_Mosmann)

clus_FlowSOM_Levine_32 <- data_FlowSOM_Levine_32[, "label"]
clus_FlowSOM_Levine_13 <- data_FlowSOM_Levine_13[, "label"]
clus_FlowSOM_Nilsson <- data_FlowSOM_Nilsson[, "label"]
clus_FlowSOM_Mosmann <- data_FlowSOM_Mosmann[, "label"]

length(clus_FlowSOM_Levine_32)
length(clus_FlowSOM_Levine_13)
length(clus_FlowSOM_Nilsson)
length(clus_FlowSOM_Mosmann)

# contingency tables

table(clus_FlowSOM_Levine_32, clus_truth_Levine_32)
table(clus_FlowSOM_Levine_13, clus_truth_Levine_13)
table(clus_FlowSOM_Nilsson, clus_truth_Nilsson)
table(clus_FlowSOM_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_FlowSOM_Levine_32 <- table(clus_FlowSOM_Levine_32)
tbl_FlowSOM_Levine_13 <- table(clus_FlowSOM_Levine_13)
tbl_FlowSOM_Nilsson <- table(clus_FlowSOM_Nilsson)
tbl_FlowSOM_Mosmann <- table(clus_FlowSOM_Mosmann)

tbl_FlowSOM_Levine_32
tbl_FlowSOM_Levine_13
tbl_FlowSOM_Nilsson
tbl_FlowSOM_Mosmann

length(tbl_FlowSOM_Levine_32)
length(tbl_FlowSOM_Levine_13)
length(tbl_FlowSOM_Nilsson)
length(tbl_FlowSOM_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FlowSOM_Levine_32 <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine_32, clus_truth_Levine_32)
res_FlowSOM_Levine_13 <- helper_match_clusters_and_evaluate(clus_FlowSOM_Levine_13, clus_truth_Levine_13)
res_FlowSOM_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Nilsson, clus_truth_Nilsson)
res_FlowSOM_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_Mosmann, clus_truth_Mosmann)

res_FlowSOM_Levine_32
res_FlowSOM_Levine_13
res_FlowSOM_Nilsson
res_FlowSOM_Mosmann




#################################
### load FlowSOM_meta results ###
#################################

# load cluster labels

file_FlowSOM_meta_Levine_32 <- "../results/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_32.txt"
file_FlowSOM_meta_Levine_13 <- "../results/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_13.txt"
file_FlowSOM_meta_Nilsson <- "../results/FlowSOM_meta/FlowSOM_meta_labels_Nilsson_2013_HSC.txt"
file_FlowSOM_meta_Mosmann <- "../results/FlowSOM_meta/FlowSOM_meta_labels_Mosmann_2014_activ.txt"

data_FlowSOM_meta_Levine_32 <- read.table(file_FlowSOM_meta_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_meta_Levine_13 <- read.table(file_FlowSOM_meta_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_meta_Nilsson <- read.table(file_FlowSOM_meta_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_FlowSOM_meta_Mosmann <- read.table(file_FlowSOM_meta_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_FlowSOM_meta_Levine_32)
head(data_FlowSOM_meta_Levine_13)
head(data_FlowSOM_meta_Nilsson)
head(data_FlowSOM_meta_Mosmann)

clus_FlowSOM_meta_Levine_32 <- data_FlowSOM_meta_Levine_32[, "label"]
clus_FlowSOM_meta_Levine_13 <- data_FlowSOM_meta_Levine_13[, "label"]
clus_FlowSOM_meta_Nilsson <- data_FlowSOM_meta_Nilsson[, "label"]
clus_FlowSOM_meta_Mosmann <- data_FlowSOM_meta_Mosmann[, "label"]

length(clus_FlowSOM_meta_Levine_32)
length(clus_FlowSOM_meta_Levine_13)
length(clus_FlowSOM_meta_Nilsson)
length(clus_FlowSOM_meta_Mosmann)

# contingency tables

table(clus_FlowSOM_meta_Levine_32, clus_truth_Levine_32)
table(clus_FlowSOM_meta_Levine_13, clus_truth_Levine_13)
table(clus_FlowSOM_meta_Nilsson, clus_truth_Nilsson)
table(clus_FlowSOM_meta_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_FlowSOM_meta_Levine_32 <- table(clus_FlowSOM_meta_Levine_32)
tbl_FlowSOM_meta_Levine_13 <- table(clus_FlowSOM_meta_Levine_13)
tbl_FlowSOM_meta_Nilsson <- table(clus_FlowSOM_meta_Nilsson)
tbl_FlowSOM_meta_Mosmann <- table(clus_FlowSOM_meta_Mosmann)

tbl_FlowSOM_meta_Levine_32
tbl_FlowSOM_meta_Levine_13
tbl_FlowSOM_meta_Nilsson
tbl_FlowSOM_meta_Mosmann

length(tbl_FlowSOM_meta_Levine_32)
length(tbl_FlowSOM_meta_Levine_13)
length(tbl_FlowSOM_meta_Nilsson)
length(tbl_FlowSOM_meta_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_FlowSOM_meta_Levine_32 <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine_32, clus_truth_Levine_32)
res_FlowSOM_meta_Levine_13 <- helper_match_clusters_and_evaluate(clus_FlowSOM_meta_Levine_13, clus_truth_Levine_13)
res_FlowSOM_meta_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Nilsson, clus_truth_Nilsson)
res_FlowSOM_meta_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_FlowSOM_meta_Mosmann, clus_truth_Mosmann)

res_FlowSOM_meta_Levine_32
res_FlowSOM_meta_Levine_13
res_FlowSOM_meta_Nilsson
res_FlowSOM_meta_Mosmann




################################
### load immunoClust results ###
################################

# load cluster labels

file_immunoClust_Levine_32 <- "../results/immunoClust/immunoClust_labels_Levine_2015_marrow_32.txt"
file_immunoClust_Levine_13 <- "../results/immunoClust/immunoClust_labels_Levine_2015_marrow_13.txt"
file_immunoClust_Nilsson <- "../results/immunoClust/immunoClust_labels_Nilsson_2013_HSC.txt"
file_immunoClust_Mosmann <- "../results/immunoClust/immunoClust_labels_Mosmann_2014_activ.txt"

data_immunoClust_Levine_32 <- read.table(file_immunoClust_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_Levine_13 <- read.table(file_immunoClust_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_Nilsson <- read.table(file_immunoClust_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_Mosmann <- read.table(file_immunoClust_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_immunoClust_Levine_32)
head(data_immunoClust_Levine_13)
head(data_immunoClust_Nilsson)
head(data_immunoClust_Mosmann)

clus_immunoClust_Levine_32 <- data_immunoClust_Levine_32[, "label"]
clus_immunoClust_Levine_13 <- data_immunoClust_Levine_13[, "label"]
clus_immunoClust_Nilsson <- data_immunoClust_Nilsson[, "label"]
clus_immunoClust_Mosmann <- data_immunoClust_Mosmann[, "label"]

length(clus_immunoClust_Levine_32)
length(clus_immunoClust_Levine_13)
length(clus_immunoClust_Nilsson)
length(clus_immunoClust_Mosmann)

# contingency tables

table(clus_immunoClust_Levine_32, clus_truth_Levine_32)
table(clus_immunoClust_Levine_13, clus_truth_Levine_13)
table(clus_immunoClust_Nilsson, clus_truth_Nilsson)
table(clus_immunoClust_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_immunoClust_Levine_32 <- table(clus_immunoClust_Levine_32)
tbl_immunoClust_Levine_13 <- table(clus_immunoClust_Levine_13)
tbl_immunoClust_Nilsson <- table(clus_immunoClust_Nilsson)
tbl_immunoClust_Mosmann <- table(clus_immunoClust_Mosmann)

tbl_immunoClust_Levine_32
tbl_immunoClust_Levine_13
tbl_immunoClust_Nilsson
tbl_immunoClust_Mosmann

length(tbl_immunoClust_Levine_32)
length(tbl_immunoClust_Levine_13)
length(tbl_immunoClust_Nilsson)
length(tbl_immunoClust_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_immunoClust_Levine_32 <- helper_match_clusters_and_evaluate(clus_immunoClust_Levine_32, clus_truth_Levine_32)
res_immunoClust_Levine_13 <- helper_match_clusters_and_evaluate(clus_immunoClust_Levine_13, clus_truth_Levine_13)
res_immunoClust_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_Nilsson, clus_truth_Nilsson)
res_immunoClust_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_Mosmann, clus_truth_Mosmann)

res_immunoClust_Levine_32
res_immunoClust_Levine_13
res_immunoClust_Nilsson
res_immunoClust_Mosmann




####################################
### load immunoClust_all results ###
####################################

# load cluster labels

file_immunoClust_all_Levine_32 <- "../results/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_32.txt"
file_immunoClust_all_Levine_13 <- "../results/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_13.txt"
file_immunoClust_all_Nilsson <- "../results/immunoClust_all/immunoClust_all_labels_Nilsson_2013_HSC.txt"
file_immunoClust_all_Mosmann <- "../results/immunoClust_all/immunoClust_all_labels_Mosmann_2014_activ.txt"

data_immunoClust_all_Levine_32 <- read.table(file_immunoClust_all_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_all_Levine_13 <- read.table(file_immunoClust_all_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_all_Nilsson <- read.table(file_immunoClust_all_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_immunoClust_all_Mosmann <- read.table(file_immunoClust_all_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_immunoClust_all_Levine_32)
head(data_immunoClust_all_Levine_13)
head(data_immunoClust_all_Nilsson)
head(data_immunoClust_all_Mosmann)

clus_immunoClust_all_Levine_32 <- data_immunoClust_all_Levine_32[, "label"]
clus_immunoClust_all_Levine_13 <- data_immunoClust_all_Levine_13[, "label"]
clus_immunoClust_all_Nilsson <- data_immunoClust_all_Nilsson[, "label"]
clus_immunoClust_all_Mosmann <- data_immunoClust_all_Mosmann[, "label"]

length(clus_immunoClust_all_Levine_32)
length(clus_immunoClust_all_Levine_13)
length(clus_immunoClust_all_Nilsson)
length(clus_immunoClust_all_Mosmann)

# contingency tables

table(clus_immunoClust_all_Levine_32, clus_truth_Levine_32)
table(clus_immunoClust_all_Levine_13, clus_truth_Levine_13)
table(clus_immunoClust_all_Nilsson, clus_truth_Nilsson)
table(clus_immunoClust_all_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_immunoClust_all_Levine_32 <- table(clus_immunoClust_all_Levine_32)
tbl_immunoClust_all_Levine_13 <- table(clus_immunoClust_all_Levine_13)
tbl_immunoClust_all_Nilsson <- table(clus_immunoClust_all_Nilsson)
tbl_immunoClust_all_Mosmann <- table(clus_immunoClust_all_Mosmann)

tbl_immunoClust_all_Levine_32
tbl_immunoClust_all_Levine_13
tbl_immunoClust_all_Nilsson
tbl_immunoClust_all_Mosmann

length(tbl_immunoClust_all_Levine_32)
length(tbl_immunoClust_all_Levine_13)
length(tbl_immunoClust_all_Nilsson)
length(tbl_immunoClust_all_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_immunoClust_all_Levine_32 <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine_32, clus_truth_Levine_32)
res_immunoClust_all_Levine_13 <- helper_match_clusters_and_evaluate(clus_immunoClust_all_Levine_13, clus_truth_Levine_13)
res_immunoClust_all_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Nilsson, clus_truth_Nilsson)
res_immunoClust_all_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_immunoClust_all_Mosmann, clus_truth_Mosmann)

res_immunoClust_all_Levine_32
res_immunoClust_all_Levine_13
res_immunoClust_all_Nilsson
res_immunoClust_all_Mosmann




############################
### load k-means results ###
############################

# load cluster labels

file_kmeans_Levine_32 <- "../results/kmeans/kmeans_labels_Levine_2015_marrow_32.txt"
file_kmeans_Levine_13 <- "../results/kmeans/kmeans_labels_Levine_2015_marrow_13.txt"
file_kmeans_Nilsson <- "../results/kmeans/kmeans_labels_Nilsson_2013_HSC.txt"
file_kmeans_Mosmann <- "../results/kmeans/kmeans_labels_Mosmann_2014_activ.txt"

data_kmeans_Levine_32 <- read.table(file_kmeans_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_kmeans_Levine_13 <- read.table(file_kmeans_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_kmeans_Nilsson <- read.table(file_kmeans_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_kmeans_Mosmann <- read.table(file_kmeans_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_kmeans_Levine_32)
head(data_kmeans_Levine_13)
head(data_kmeans_Nilsson)
head(data_kmeans_Mosmann)

clus_kmeans_Levine_32 <- data_kmeans_Levine_32[, "label"]
clus_kmeans_Levine_13 <- data_kmeans_Levine_13[, "label"]
clus_kmeans_Nilsson <- data_kmeans_Nilsson[, "label"]
clus_kmeans_Mosmann <- data_kmeans_Mosmann[, "label"]

length(clus_kmeans_Levine_32)
length(clus_kmeans_Levine_13)
length(clus_kmeans_Nilsson)
length(clus_kmeans_Mosmann)

# contingency tables

table(clus_kmeans_Levine_32, clus_truth_Levine_32)
table(clus_kmeans_Levine_13, clus_truth_Levine_13)
table(clus_kmeans_Nilsson, clus_truth_Nilsson)
table(clus_kmeans_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_kmeans_Levine_32 <- table(clus_kmeans_Levine_32)
tbl_kmeans_Levine_13 <- table(clus_kmeans_Levine_13)
tbl_kmeans_Nilsson <- table(clus_kmeans_Nilsson)
tbl_kmeans_Mosmann <- table(clus_kmeans_Mosmann)

tbl_kmeans_Levine_32
tbl_kmeans_Levine_13
tbl_kmeans_Nilsson
tbl_kmeans_Mosmann

length(tbl_kmeans_Levine_32)
length(tbl_kmeans_Levine_13)
length(tbl_kmeans_Nilsson)
length(tbl_kmeans_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_kmeans_Levine_32 <- helper_match_clusters_and_evaluate(clus_kmeans_Levine_32, clus_truth_Levine_32)
res_kmeans_Levine_13 <- helper_match_clusters_and_evaluate(clus_kmeans_Levine_13, clus_truth_Levine_13)
res_kmeans_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_kmeans_Nilsson, clus_truth_Nilsson)
res_kmeans_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_kmeans_Mosmann, clus_truth_Mosmann)

res_kmeans_Levine_32
res_kmeans_Levine_13
res_kmeans_Nilsson
res_kmeans_Mosmann




###############################
### load Rclusterpp results ###
###############################

# load cluster labels

file_Rclusterpp_Levine_32 <- "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_32.txt"
file_Rclusterpp_Levine_13 <- "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_13.txt"
file_Rclusterpp_Nilsson <- "../results/Rclusterpp/Rclusterpp_labels_Nilsson_2013_HSC.txt"
# file_Rclusterpp_Mosmann <- "../results/Rclusterpp/Rclusterpp_labels_Mosmann_2014_activ.txt"

data_Rclusterpp_Levine_32 <- read.table(file_Rclusterpp_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_Rclusterpp_Levine_13 <- read.table(file_Rclusterpp_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_Rclusterpp_Nilsson <- read.table(file_Rclusterpp_Nilsson, header = TRUE, sep = "\t", comment.char = "")
# data_Rclusterpp_Mosmann <- read.table(file_Rclusterpp_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_Rclusterpp_Levine_32)
head(data_Rclusterpp_Levine_13)
head(data_Rclusterpp_Nilsson)
# head(data_Rclusterpp_Mosmann)

clus_Rclusterpp_Levine_32 <- data_Rclusterpp_Levine_32[, "label"]
clus_Rclusterpp_Levine_13 <- data_Rclusterpp_Levine_13[, "label"]
clus_Rclusterpp_Nilsson <- data_Rclusterpp_Nilsson[, "label"]
# clus_Rclusterpp_Mosmann <- data_Rclusterpp_Mosmann[, "label"]

length(clus_Rclusterpp_Levine_32)
length(clus_Rclusterpp_Levine_13)
length(clus_Rclusterpp_Nilsson)
# length(clus_Rclusterpp_Mosmann)

# contingency tables

table(clus_Rclusterpp_Levine_32, clus_truth_Levine_32)
table(clus_Rclusterpp_Levine_13, clus_truth_Levine_13)
table(clus_Rclusterpp_Nilsson, clus_truth_Nilsson)
# table(clus_Rclusterpp_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_Rclusterpp_Levine_32 <- table(clus_Rclusterpp_Levine_32)
tbl_Rclusterpp_Levine_13 <- table(clus_Rclusterpp_Levine_13)
tbl_Rclusterpp_Nilsson <- table(clus_Rclusterpp_Nilsson)
# tbl_Rclusterpp_Mosmann <- table(clus_Rclusterpp_Mosmann)

tbl_Rclusterpp_Levine_32
tbl_Rclusterpp_Levine_13
tbl_Rclusterpp_Nilsson
# tbl_Rclusterpp_Mosmann

length(tbl_Rclusterpp_Levine_32)
length(tbl_Rclusterpp_Levine_13)
length(tbl_Rclusterpp_Nilsson)
# length(tbl_Rclusterpp_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_Rclusterpp_Levine_32 <- helper_match_clusters_and_evaluate(clus_Rclusterpp_Levine_32, clus_truth_Levine_32)
res_Rclusterpp_Levine_13 <- helper_match_clusters_and_evaluate(clus_Rclusterpp_Levine_13, clus_truth_Levine_13)
res_Rclusterpp_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_Rclusterpp_Nilsson, clus_truth_Nilsson)
# res_Rclusterpp_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_Rclusterpp_Mosmann, clus_truth_Mosmann)

res_Rclusterpp_Levine_32
res_Rclusterpp_Levine_13
res_Rclusterpp_Nilsson
# res_Rclusterpp_Mosmann




################################
### load SamSPECTRAL results ###
################################

# load cluster labels

file_SamSPECTRAL_Levine_32 <- "../results/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_32.txt"
file_SamSPECTRAL_Levine_13 <- "../results/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_13.txt"
file_SamSPECTRAL_Nilsson <- "../results/SamSPECTRAL/SamSPECTRAL_labels_Nilsson_2013_HSC.txt"
file_SamSPECTRAL_Mosmann <- "../results/SamSPECTRAL/SamSPECTRAL_labels_Mosmann_2014_activ.txt"

data_SamSPECTRAL_Levine_32 <- read.table(file_SamSPECTRAL_Levine_32, header = TRUE, sep = "\t", comment.char = "")
data_SamSPECTRAL_Levine_13 <- read.table(file_SamSPECTRAL_Levine_13, header = TRUE, sep = "\t", comment.char = "")
data_SamSPECTRAL_Nilsson <- read.table(file_SamSPECTRAL_Nilsson, header = TRUE, sep = "\t", comment.char = "")
data_SamSPECTRAL_Mosmann <- read.table(file_SamSPECTRAL_Mosmann, header = TRUE, sep = "\t", comment.char = "")

head(data_SamSPECTRAL_Levine_32)
head(data_SamSPECTRAL_Levine_13)
head(data_SamSPECTRAL_Nilsson)
head(data_SamSPECTRAL_Mosmann)

clus_SamSPECTRAL_Levine_32 <- data_SamSPECTRAL_Levine_32[, "label"]
clus_SamSPECTRAL_Levine_13 <- data_SamSPECTRAL_Levine_13[, "label"]
clus_SamSPECTRAL_Nilsson <- data_SamSPECTRAL_Nilsson[, "label"]
clus_SamSPECTRAL_Mosmann <- data_SamSPECTRAL_Mosmann[, "label"]

length(clus_SamSPECTRAL_Levine_32)
length(clus_SamSPECTRAL_Levine_13)
length(clus_SamSPECTRAL_Nilsson)
length(clus_SamSPECTRAL_Mosmann)

# contingency tables

table(clus_SamSPECTRAL_Levine_32, clus_truth_Levine_32)
table(clus_SamSPECTRAL_Levine_13, clus_truth_Levine_13)
table(clus_SamSPECTRAL_Nilsson, clus_truth_Nilsson)
table(clus_SamSPECTRAL_Mosmann, clus_truth_Mosmann)

# cluster sizes and number of clusters

tbl_SamSPECTRAL_Levine_32 <- table(clus_SamSPECTRAL_Levine_32)
tbl_SamSPECTRAL_Levine_13 <- table(clus_SamSPECTRAL_Levine_13)
tbl_SamSPECTRAL_Nilsson <- table(clus_SamSPECTRAL_Nilsson)
tbl_SamSPECTRAL_Mosmann <- table(clus_SamSPECTRAL_Mosmann)

tbl_SamSPECTRAL_Levine_32
tbl_SamSPECTRAL_Levine_13
tbl_SamSPECTRAL_Nilsson
tbl_SamSPECTRAL_Mosmann

length(tbl_SamSPECTRAL_Levine_32)
length(tbl_SamSPECTRAL_Levine_13)
length(tbl_SamSPECTRAL_Nilsson)
length(tbl_SamSPECTRAL_Mosmann)

# match cluster labels by highest F1 score and calculate results
# precision, recall, F1 score, matched cluster labels, number of cells per matched cluster

res_SamSPECTRAL_Levine_32 <- helper_match_clusters_and_evaluate(clus_SamSPECTRAL_Levine_32, clus_truth_Levine_32)
res_SamSPECTRAL_Levine_13 <- helper_match_clusters_and_evaluate(clus_SamSPECTRAL_Levine_13, clus_truth_Levine_13)
res_SamSPECTRAL_Nilsson <- helper_match_one_rare_cluster_and_evaluate(clus_SamSPECTRAL_Nilsson, clus_truth_Nilsson)
res_SamSPECTRAL_Mosmann <- helper_match_one_rare_cluster_and_evaluate(clus_SamSPECTRAL_Mosmann, clus_truth_Mosmann)

res_SamSPECTRAL_Levine_32
res_SamSPECTRAL_Levine_13
res_SamSPECTRAL_Nilsson
res_SamSPECTRAL_Mosmann


