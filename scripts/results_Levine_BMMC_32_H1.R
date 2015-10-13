# R script to compare clustering results from different methods, for the data set
# "Levine_BMMC_32_H1".
#
# Lukas M. Weber, October 2015


library(flowCore)


# load manually gated cluster labels ("truth")

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

file_PhenoGraph <- "../results/Levine_BMMC_32_H1/Levine_BMMC_32_H1_PhenoGraph.fcs"
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


