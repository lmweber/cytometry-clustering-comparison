#########################################################################################
# R script to load results for each clustering method
#
# data set "Levine_BMMC_32"
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)

source("match_clusters_and_evaluate.R")


# =================================
# manually gated clusters ("truth")
# =================================

file_truth <- "../data/Levine_BMMC_32/Levine_BMMC_32.fcs"
data_truth <- flowCore::exprs(flowCore::read.FCS(file_truth, transformation = FALSE))

head(data_truth)
dim(data_truth)

clus_truth <- data_truth[, "label"]
length(clus_truth)

tbl_truth <- table(clus_truth)  # cluster sizes
tbl_truth
n_clus_truth <- length(tbl_truth)  # number of clusters
n_clus_truth

n_cells_truth <- tbl_truth  # number of cells per cluster



# ===============
# ACCENSE results
# ===============

# load cluster labels

file_ACCENSE <- "../results/Levine_BMMC_32/ACCENSE/accense_output.csv"
data_ACCENSE <- read.csv(file_ACCENSE, header = TRUE)

head(data_ACCENSE)
dim(data_ACCENSE)  # note this is a subsample of the full data set

# true cluster labels for subsampled points

clus_truth_ACCENSE <- data_ACCENSE[, "label"]
table(clus_truth_ACCENSE)

# ACCENSE cluster labels

clus_ACCENSE <- data_ACCENSE[, "population"]

table(clus_ACCENSE, clus_truth_ACCENSE)  # contingency table

tbl_ACCENSE <- table(clus_ACCENSE)  # cluster sizes
tbl_ACCENSE
n_clus_ACCENSE <- length(tbl_ACCENSE)  # number of clusters
n_clus_ACCENSE

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_ACCENSE <- match_clusters_and_evaluate(clus_ACCENSE, clus_truth_ACCENSE)

pr_ACCENSE <- res_ACCENSE$pr
re_ACCENSE <- res_ACCENSE$re
F1_ACCENSE <- res_ACCENSE$F1

# matched cluster labels

labels_matched_ACCENSE <- res_ACCENSE$labels_matched

# number of cells per matched cluster

n_cells_ACCENSE <- res_ACCENSE$n_cells



# ==============
# DensVM results
# ==============

# load cluster labels

file_DensVM <- "../results/Levine_BMMC_32/DensVM/cytofkit_analysis_tsne_cluster.txt"
data_DensVM <- read.table(file_DensVM, header = TRUE, sep = "\t")

head(data_DensVM)
dim(data_DensVM)  # note this is a subsample of the full data set

# true cluster labels for subsampled points

file_truth_DensVM <- "../results/Levine_BMMC_32/DensVM/cytofkit_analysis_analyzedFCS/Levine_BMMC_32_notransf.fcs"
data_truth_DensVM <- flowCore::exprs(flowCore::read.FCS(file_truth_DensVM, transformation = FALSE))

head(data_truth_DensVM)
dim(data_truth_DensVM)

clus_truth_DensVM <- data_truth_DensVM[, "label"]

table(clus_truth_DensVM)
length(clus_truth_DensVM)
length(table(clus_truth_DensVM))

# DensVM cluster labels

clus_DensVM <- data_DensVM[, "cluster"]

table(clus_DensVM, clus_truth_DensVM)  # contingency table

tbl_DensVM <- table(clus_DensVM)  # cluster sizes
tbl_DensVM
n_clus_DensVM <- length(tbl_DensVM)  # number of clusters
n_clus_DensVM

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_DensVM <- match_clusters_and_evaluate(clus_DensVM, clus_truth_DensVM)

pr_DensVM <- res_DensVM$pr
re_DensVM <- res_DensVM$re
F1_DensVM <- res_DensVM$F1

# matched cluster labels

labels_matched_DensVM <- res_DensVM$labels_matched

# number of cells per matched cluster

n_cells_DensVM <- res_DensVM$n_cells



# ===============
# FlowSOM results
# ===============

# load cluster labels

file_FlowSOM <- "../results/Levine_BMMC_32/FlowSOM/Levine_BMMC_32_FlowSOM_labels.txt"
data_FlowSOM <- read.table(file_FlowSOM, header = TRUE, sep = "\t", comment.char = "")

head(data_FlowSOM)
dim(data_FlowSOM)

clus_FlowSOM <- data_FlowSOM[, "label"]

table(clus_FlowSOM, clus_truth)  # contingency table

tbl_FlowSOM <- table(clus_FlowSOM)  # cluster sizes
tbl_FlowSOM
n_clus_FlowSOM <- length(tbl_FlowSOM)  # number of clusters
n_clus_FlowSOM

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_FlowSOM <- match_clusters_and_evaluate(clus_FlowSOM, clus_truth)

pr_FlowSOM <- res_FlowSOM$pr
re_FlowSOM <- res_FlowSOM$re
F1_FlowSOM <- res_FlowSOM$F1

# matched cluster labels

labels_matched_FlowSOM <- res_FlowSOM$labels_matched

# number of cells per matched cluster

n_cells_FlowSOM <- res_FlowSOM$n_cells



# ====================
# FlowSOM_meta results
# ====================

# load cluster labels

file_FlowSOM_meta <- "../results/Levine_BMMC_32/FlowSOM_meta/Levine_BMMC_32_FlowSOM_meta_labels.txt"
data_FlowSOM_meta <- read.table(file_FlowSOM_meta, header = TRUE, sep = "\t", comment.char = "")

head(data_FlowSOM_meta)
dim(data_FlowSOM_meta)

clus_FlowSOM_meta <- data_FlowSOM_meta[, "label"]

table(clus_FlowSOM_meta, clus_truth)  # contingency table

tbl_FlowSOM_meta <- table(clus_FlowSOM_meta)  # cluster sizes
tbl_FlowSOM_meta
n_clus_FlowSOM_meta <- length(tbl_FlowSOM_meta)  # number of clusters
n_clus_FlowSOM_meta

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_FlowSOM_meta <- match_clusters_and_evaluate(clus_FlowSOM_meta, clus_truth)

pr_FlowSOM_meta <- res_FlowSOM_meta$pr
re_FlowSOM_meta <- res_FlowSOM_meta$re
F1_FlowSOM_meta <- res_FlowSOM_meta$F1

# matched cluster labels

labels_matched_FlowSOM_meta <- res_FlowSOM_meta$labels_matched

# number of cells per matched cluster

n_cells_FlowSOM_meta <- res_FlowSOM_meta$n_cells



# ===================
# immunoClust results
# ===================

# load cluster labels

file_immunoClust <- "../results/Levine_BMMC_32/immunoClust/Levine_BMMC_32_immunoClust_labels.txt"
data_immunoClust <- read.table(file_immunoClust, header = TRUE, sep = "\t")

head(data_immunoClust)
dim(data_immunoClust)

clus_immunoClust <- data_immunoClust[, "label"]

table(clus_immunoClust, clus_truth)  # contingency table

tbl_immunoClust <- table(clus_immunoClust)  # cluster sizes
tbl_immunoClust
n_clus_immunoClust <- length(tbl_immunoClust)  # number of clusters
n_clus_immunoClust

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_immunoClust <- match_clusters_and_evaluate(clus_immunoClust, clus_truth)

pr_immunoClust <- res_immunoClust$pr
re_immunoClust <- res_immunoClust$re
F1_immunoClust <- res_immunoClust$F1

# matched cluster labels

labels_matched_immunoClust <- res_immunoClust$labels_matched

# number of cells per matched cluster

n_cells_immunoClust <- res_immunoClust$n_cells



# =======================
# immunoClust_all results
# =======================

# load cluster labels

file_immunoClust_all <- "../results/Levine_BMMC_32/immunoClust_all/Levine_BMMC_32_immunoClust_all_labels.txt"
data_immunoClust_all <- read.table(file_immunoClust_all, header = TRUE, sep = "\t")

head(data_immunoClust_all)
dim(data_immunoClust_all)

clus_immunoClust_all <- data_immunoClust_all[, "label"]

table(clus_immunoClust_all, clus_truth)  # contingency table

tbl_immunoClust_all <- table(clus_immunoClust_all)  # cluster sizes
tbl_immunoClust_all
n_clus_immunoClust_all <- length(tbl_immunoClust_all)  # number of cluster
n_clus_immunoClust_all

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_immunoClust_all <- match_clusters_and_evaluate(clus_immunoClust_all, clus_truth)

pr_immunoClust_all <- res_immunoClust_all$pr
re_immunoClust_all <- res_immunoClust_all$re
F1_immunoClust_all <- res_immunoClust_all$F1

# matched cluster labels

labels_matched_immunoClust_all <- res_immunoClust_all$labels_matched

# number of cells per matched cluster

n_cells_immunoClust_all <- res_immunoClust_all$n_cells



# ===============
# k-means results
# ===============

# load cluster labels

file_kmeans <- "../results/Levine_BMMC_32/kmeans/Levine_BMMC_32_kmeans_labels.txt"
data_kmeans <- read.table(file_kmeans, header = TRUE, sep = "\t")

head(data_kmeans)
dim(data_kmeans)

clus_kmeans <- data_kmeans[, "label"]

table(clus_kmeans, clus_truth)  # contingency table

tbl_kmeans <- table(clus_kmeans)  # cluster sizes
tbl_kmeans
n_clus_kmeans <- length(tbl_kmeans)  # number of clusters
n_clus_kmeans

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_kmeans <- match_clusters_and_evaluate(clus_kmeans, clus_truth)

pr_kmeans <- res_kmeans$pr
re_kmeans <- res_kmeans$re
F1_kmeans <- res_kmeans$F1

# matched cluster labels

labels_matched_kmeans <- res_kmeans$labels_matched

# number of cells per matched cluster

n_cells_kmeans <- res_kmeans$n_cells



# ==================
# PhenoGraph results
# ==================

# load cluster labels

file_PhenoGraph <- "../results/Levine_BMMC_32/PhenoGraph/Levine_BMMC_32_PhenoGraph.fcs"
data_PhenoGraph <- flowCore::exprs(flowCore::read.FCS(file_PhenoGraph, transformation = FALSE))

head(data_PhenoGraph)
dim(data_PhenoGraph)

clus_PhenoGraph <- data_PhenoGraph[, grep("PhenoGraph", colnames(data_PhenoGraph))]

table(clus_PhenoGraph, clus_truth)  # contingency table

tbl_PhenoGraph <- table(clus_PhenoGraph)  # cluster sizes
tbl_PhenoGraph
n_clus_PhenoGraph <- length(tbl_PhenoGraph)  # number of clusters
n_clus_PhenoGraph

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_PhenoGraph <- match_clusters_and_evaluate(clus_PhenoGraph, clus_truth)

pr_PhenoGraph <- res_PhenoGraph$pr
re_PhenoGraph <- res_PhenoGraph$re
F1_PhenoGraph <- res_PhenoGraph$F1

# matched cluster labels

labels_matched_PhenoGraph <- res_PhenoGraph$labels_matched

# number of cells per matched cluster

n_cells_PhenoGraph <- res_PhenoGraph$n_cells



# ==================
# Rclusterpp results
# ==================

# load cluster labels

file_Rclusterpp <- "../results/Levine_BMMC_32/Rclusterpp/Levine_BMMC_32_Rclusterpp_labels.txt"
data_Rclusterpp <- read.table(file_Rclusterpp, header = TRUE, sep = "\t")

head(data_Rclusterpp)
dim(data_Rclusterpp)

clus_Rclusterpp <- data_Rclusterpp[, "label"]

table(clus_Rclusterpp, clus_truth)  # contingency table

tbl_Rclusterpp <- table(clus_Rclusterpp)  # cluster sizes
tbl_Rclusterpp
n_clus_Rclusterpp <- length(tbl_Rclusterpp)  # number of clusters
n_clus_Rclusterpp

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_Rclusterpp <- match_clusters_and_evaluate(clus_Rclusterpp, clus_truth)

pr_Rclusterpp <- res_Rclusterpp$pr
re_Rclusterpp <- res_Rclusterpp$re
F1_Rclusterpp <- res_Rclusterpp$F1

# matched cluster labels

labels_matched_Rclusterpp <- res_Rclusterpp$labels_matched

# number of cells per matched cluster

n_cells_Rclusterpp <- res_Rclusterpp$n_cells



# =============
# SWIFT results
# =============

# load cluster labels

file_SWIFT <- "../results/Levine_BMMC_32/SWIFT/Levine_BMMC_32_notransf.fcs.Cluster_Output.txt"
data_SWIFT <- read.table(file_SWIFT, header = TRUE, sep = "\t", comment.char = "")

head(data_SWIFT)
dim(data_SWIFT)

clus_SWIFT <- data_SWIFT[, "MergeCluster."]

table(clus_SWIFT, clus_truth)  # contingency table

tbl_SWIFT <- table(clus_SWIFT)  # cluster sizes
tbl_SWIFT
n_clus_SWIFT <- length(tbl_SWIFT)  # number of clusters
n_clus_SWIFT

# match cluster labels by highest F1 score; and calculate precision, recall, and F1 score

res_SWIFT <- match_clusters_and_evaluate(clus_SWIFT, clus_truth)

pr_SWIFT <- res_SWIFT$pr
re_SWIFT <- res_SWIFT$re
F1_SWIFT <- res_SWIFT$F1

# matched cluster labels

labels_matched_SWIFT <- res_SWIFT$labels_matched

# number of cells per matched cluster

n_cells_SWIFT <- res_SWIFT$n_cells


