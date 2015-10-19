# R script to compare clustering results for data set "Levine_BMMC_32_H1"
#
# Lukas Weber, October 2015


library(flowCore)


# load manually gated cluster labels ("truth")
# --------------------------------------------

file_truth <- "../data/Levine_BMMC_32_H1/Levine_BMMC_32_H1.txt"
data_truth <- read.table(file_truth, header = TRUE, sep = "\t")
head(data_truth)
dim(data_truth)

clus_truth <- data_truth$label
length(clus_truth)

max(clus_truth)  # total number of clusters (including "unassigned" group)

table(clus_truth)

res <- data.frame(truth = clus_truth)



# load PhenoGraph clusters
# ------------------------

file_PhenoGraph <- "../results/Levine_BMMC_32_H1/PhenoGraph/Levine_BMMC_32_H1_PhenoGraph.fcs"
data_PhenoGraph <- flowCore::exprs(flowCore::read.FCS(file_PhenoGraph, transformation = FALSE))
head(data_PhenoGraph)
dim(data_PhenoGraph)

clus_PhenoGraph <- data_PhenoGraph[, grep("PhenoGraph", colnames(data_PhenoGraph))]
head(clus_PhenoGraph)
length(clus_PhenoGraph)

max(clus_PhenoGraph)  # total number of clusters detected

table(clus_PhenoGraph)


# match PhenoGraph clusters to manually gated clusters
# by ordering rows of table

tbl_PhenoGraph <- table(clus_PhenoGraph, clus_truth)
tbl_PhenoGraph

# sort table
tbl_PhenoGraph_sorted <- matrix(nrow = 0, ncol = ncol(tbl_PhenoGraph))
tbl_PhenoGraph_temp <- tbl_PhenoGraph
for (i in 1:ncol(tbl_PhenoGraph_temp)) {
  ix <- which.max(tbl_PhenoGraph_temp[, i])
  tbl_PhenoGraph_sorted <- rbind(tbl_PhenoGraph_sorted, tbl_PhenoGraph_temp[ix, ])
  rownames(tbl_PhenoGraph_sorted)[i] <- rownames(tbl_PhenoGraph_temp)[ix]
  tbl_PhenoGraph_temp <- tbl_PhenoGraph_temp[-ix, ]
}
tbl_PhenoGraph_sorted <- rbind(tbl_PhenoGraph_sorted, tbl_PhenoGraph_temp)
tbl_PhenoGraph_sorted

# matched cluster labels
#clus_PhenoGraph_matched <- factor(clus_PhenoGraph, 
#                                  levels = as.integer(rownames(tbl_PhenoGraph_sorted)))
#levels(clus_PhenoGraph_matched) <- 1:nrow(tbl_PhenoGraph_sorted)

#head(clus_PhenoGraph_matched)
#length(clus_PhenoGraph_matched)
#table(clus_PhenoGraph_matched, clus_truth)

# store results
#res <- cbind(res, PhenoGraph = clus_PhenoGraph_matched)


# calculate F1-measure for each cluster

F1_PhenoGraph <- rep(NA, ncol(tbl_PhenoGraph_sorted) - 1)
for (i in seq_along(F1_PhenoGraph)) {
  precision_i <- tbl_PhenoGraph_sorted[i, i] / sum(tbl_PhenoGraph_sorted[i, -ncol(tbl_PhenoGraph_sorted)])  # drop unassigned
  recall_i <- tbl_PhenoGraph_sorted[i, i] / sum(tbl_PhenoGraph_sorted[, i])
  F1_PhenoGraph[i] <- 2 * (precision_i * recall_i) / (precision_i + recall_i)
}
F1_PhenoGraph
F1_PhenoGraph[F1_PhenoGraph == "NaN"] <- 0  # replace with zeros
F1_PhenoGraph

# match cluster names

names(F1_PhenoGraph) <- colnames(tbl_PhenoGraph_sorted)[-ncol(tbl_PhenoGraph_sorted)]
F1_PhenoGraph



# load DensVM clusters
# --------------------

file_DensVM <- "../results/Levine_BMMC_32_H1/DensVM/cytof_tsne_cluster.txt"
data_DensVM <- read.table(file_DensVM, header = TRUE, sep = "\t")
head(data_DensVM)
dim(data_DensVM)

clus_DensVM <- data_DensVM$cluster
length(clus_DensVM)
table(clus_DensVM)
max(clus_DensVM)

# load truth labels for the 20,000 subsampled points

file_truth_DensVM <- "../results/Levine_BMMC_32_H1/DensVM/cytof_analyzedFCS/Levine_BMMC_32_H1_notransform.fcs"
data_truth_DensVM <- flowCore::exprs(flowCore::read.FCS(file_truth_DensVM, transformation = FALSE))
head(data_truth_DensVM)

clus_truth_DensVM <- data_truth_DensVM[, "label"]

# compare clustering

tbl_DensVM <- table(clus_DensVM, clus_truth_DensVM)  # note fewer clusters, i.e. DensVM has lower resolution than manual gating
tbl_DensVM

# rearrange to match clusters

# sort table
tbl_DensVM_sorted <- matrix(nrow = nrow(tbl_DensVM), ncol = 0)
tbl_DensVM_temp <- tbl_DensVM
for (i in 1:nrow(tbl_DensVM_temp)) {  # sort by rows since nrow < ncol
  ix <- which.max(tbl_DensVM_temp[i, -ncol(tbl_DensVM_temp)])  # drop unassigned column
  tbl_DensVM_sorted <- cbind(tbl_DensVM_sorted, tbl_DensVM_temp[, ix])
  colnames(tbl_DensVM_sorted)[i] <- colnames(tbl_DensVM_temp)[ix]
  tbl_DensVM_temp <- tbl_DensVM_temp[, -ix]
}
tbl_DensVM_sorted <- cbind(tbl_DensVM_sorted, tbl_DensVM_temp)
tbl_DensVM_sorted

# calculate F1-measure for each cluster

F1_DensVM <- rep(NA, nrow(tbl_DensVM_sorted))
for (i in seq_along(F1_DensVM)) {
  precision_i <- tbl_DensVM_sorted[i, i] / sum(tbl_DensVM_sorted[i, -ncol(tbl_DensVM_sorted)])  # drop unassigned
  recall_i <- tbl_DensVM_sorted[i, i] / sum(tbl_DensVM_sorted[, i])
  F1_DensVM[i] <- 2 * (precision_i * recall_i) / (precision_i + recall_i)
}
F1_DensVM

# match cluster names

names(F1_DensVM) <- colnames(tbl_DensVM_sorted)[1:length(F1_DensVM)]
F1_DensVM

F1_DensVM_updated <- rep(NA, ncol(tbl_DensVM_sorted) - 1)
F1_DensVM_updated[as.integer(names(F1_DensVM))] <- F1_DensVM
F1_DensVM_updated
F1_DensVM_updated[is.na(F1_DensVM_updated)] <- 0

F1_DensVM <- F1_DensVM_updated
F1_DensVM



# ACCENSE
# -------

file_ACCENSE <- "../results/Levine_BMMC_32_H1/ACCENSE/accense_output.csv"
data_ACCENSE <- read.csv(file_ACCENSE, header = TRUE)
head(data_ACCENSE)
dim(data_ACCENSE)  # note these are not the same 20,000 points as in the DensVM data frame

clus_ACCENSE <- data_ACCENSE$population
table(clus_ACCENSE)
max(clus_ACCENSE)

clus_truth_ACCENSE <- data_ACCENSE$label

tbl_ACCENSE <- table(clus_ACCENSE, clus_truth_ACCENSE)
tbl_ACCENSE

# rearrange to match clusters

# sort table (by column)
tbl_ACCENSE_sorted <- matrix(nrow = 0, ncol = ncol(tbl_ACCENSE))
tbl_ACCENSE_temp <- tbl_ACCENSE
for (i in 1:ncol(tbl_ACCENSE_temp)) {
  ix <- which.max(tbl_ACCENSE_temp[, i])
  tbl_ACCENSE_sorted <- rbind(tbl_ACCENSE_sorted, tbl_ACCENSE_temp[ix, ])
  rownames(tbl_ACCENSE_sorted)[i] <- rownames(tbl_ACCENSE_temp)[ix]
  tbl_ACCENSE_temp <- tbl_ACCENSE_temp[-ix, ]
}
tbl_ACCENSE_sorted <- rbind(tbl_ACCENSE_sorted, tbl_ACCENSE_temp)
tbl_ACCENSE_sorted

# calculate F1-measure for each cluster

F1_ACCENSE <- rep(NA, ncol(tbl_ACCENSE_sorted) - 1)
for (i in seq_along(F1_ACCENSE)) {
  precision_i <- tbl_ACCENSE_sorted[i, i] / sum(tbl_ACCENSE_sorted[i, -ncol(tbl_ACCENSE_sorted)])  # drop unassigned
  recall_i <- tbl_ACCENSE_sorted[i, i] / sum(tbl_ACCENSE_sorted[, i])
  F1_ACCENSE[i] <- 2 * (precision_i * recall_i) / (precision_i + recall_i)
}
F1_ACCENSE
F1_ACCENSE[F1_ACCENSE == "NaN"] <- 0  # replace with zeros
F1_ACCENSE

# match cluster names

names(F1_ACCENSE) <- colnames(tbl_ACCENSE_sorted)[-ncol(tbl_ACCENSE_sorted)]
F1_ACCENSE



# Plot heatmap of F1-measures for individual clusters
# ---------------------------------------------------

library(pheatmap)
library(RColorBrewer)

data_heatmap <- cbind(F1_PhenoGraph, F1_ACCENSE, F1_DensVM)
colnames(data_heatmap) <- gsub("^F1_", "", colnames(data_heatmap))

png("../plots/heatmap_F1_measure_clusters_Levine_BMMC_32_H1.png", width = 600, height = 1000)
pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, fontsize_number = 13, 
         cex = 1.5)
dev.off()



# "True" cluster sizes
# --------------------

colSums(tbl_PhenoGraph_sorted)


