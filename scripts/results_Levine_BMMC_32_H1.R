# R script to compare clustering results from different methods, for the data set
# "Levine_BMMC_32_H1".
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
# by ordering rows of confusion matrix

tbl_PhenoGraph <- table(clus_PhenoGraph, clus_truth)
tbl_PhenoGraph

# sort confusion matrix
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
clus_PhenoGraph_matched <- factor(clus_PhenoGraph, 
                                  levels = as.integer(rownames(tbl_PhenoGraph_sorted)))
levels(clus_PhenoGraph_matched) <- 1:nrow(tbl_PhenoGraph_sorted)

head(clus_PhenoGraph_matched)
length(clus_PhenoGraph_matched)
table(clus_PhenoGraph_matched, clus_truth)

# store results
res <- cbind(res, PhenoGraph = clus_PhenoGraph_matched)


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

table(clus_DensVM, clus_truth)  # require labels from subsampled data


####### alternatively

# load DensVM FCS output file

file_DensVM <- "../results/Levine_BMMC_32_H1/DensVM/cytof_analyzedFCS/Levine_BMMC_32_H1_notransform.fcs"
data_DensVM <- flowCore::exprs(flowCore::read.FCS(file_DensVM, transformation = FALSE))
head(data_DensVM)

clus_DensVM_truth <- data_DensVM[, "label"]

table(clus_DensVM, clus_DensVM_truth)
## these match very well, although some are combined (DensVM has less resolution than manual gating)



# ACCENSE
# -------

file_ACCENSE <- "../results/Levine_BMMC_32_H1/ACCENSE/accense_output.csv"
data_ACCENSE <- read.csv(file_ACCENSE, header = TRUE)
head(data_ACCENSE)

clus_ACCENSE <- data_ACCENSE$population
table(clus_ACCENSE)
max(clus_ACCENSE)

clus_ACCENSE_truth <- data_ACCENSE$label

table(clus_ACCENSE, clus_ACCENSE_truth)
## these are pretty good too
## but neither ACCENSE or DensVM will probably work well for rare cell populations

