#########################################################################################
# R script to run FlowSOM and FlowSOM_meta
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(FlowSOM)


# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32.fcs"
data <- flowCore::read.FCS(file, transformation = FALSE)  # needs to be a flowFrame object

head(data)
dim(data)

marker_cols <- (1:ncol(data))[-grep("label", colnames(data))]  # columns to use to calculate SOM


# run FlowSOM

set.seed(123)

system.time(
fSOM <- FlowSOM::ReadInput(data, transform = FALSE, scale = FALSE, silent = FALSE)
)

system.time(
fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols)  # very fast
)

system.time(
fSOM <- FlowSOM::BuildMST(fSOM)
)

FlowSOM::PlotStars(fSOM)


# extract cluster labels

str(fSOM$map)

head(fSOM$map$mapping)
dim(fSOM$map$mapping)

clus_FlowSOM <- fSOM$map$mapping[, 1]  # cluster labels

length(clus_FlowSOM)

table(clus_FlowSOM)
length(table(clus_FlowSOM))  # number of clusters


# optional metaclustering step ("FlowSOM_meta")

k <- 20  # number of clusters

system.time(
meta_clustering <- FlowSOM::metaClustering_consensus(fSOM$map$codes, k = k)
)
meta_clustering

meta_clustering_per_cell <- meta_clustering[fSOM$map$mapping[, 1]]

clus_FlowSOM_meta <- meta_clustering_per_cell

length(clus_FlowSOM_meta)

table(clus_FlowSOM_meta)
length(table(clus_FlowSOM_meta))  # number of clusters


# save files

res_FlowSOM <- data.frame(label = clus_FlowSOM)

write.table(data.frame(label = clus_FlowSOM), 
            file = "../results/Levine_BMMC_32/FlowSOM/Levine_BMMC_32_FlowSOM_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(data.frame(label = clus_FlowSOM_meta), 
            file = "../results/Levine_BMMC_32/FlowSOM_meta/Levine_BMMC_32_FlowSOM_meta_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

