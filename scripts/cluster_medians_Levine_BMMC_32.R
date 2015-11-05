#########################################################################################
# R script to calculate median marker expression values for each cluster
#
# data set "Levine_BMMC_32"
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(robustbase)
library(pheatmap)
library(RColorBrewer)


# helper functions

source("match_clusters_and_evaluate.R")
source("calculate_cluster_medians_scaled.R")


# results

source("load_results_Levine_BMMC_32.R")
source("heatmaps_tables_Levine_BMMC_32.R")


#save.image("../current.RData")
#load("../current.RData")



############################################
### Cluster medians along each dimension ###
############################################


# =======================================
# cluster medians: truth (manually gated)
# =======================================

data_medians <- as.data.frame(data_truth[, -grep("label", colnames(data_truth))])
dim(data_medians)

colnames(data_medians) <- gsub("\\(.*$", "", colnames(data_medians))

table(clus_truth)
length(clus_truth)

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_truth <- calculate_cluster_medians_scaled(data_medians, clus_truth)

# plot

heatmap_truth <- pheatmap(medians_truth, 
                          color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
                          cluster_rows = TRUE, cluster_cols = TRUE, 
                          main = "Cluster medians: manually gated", 
                          filename = "../plots/cluster_medians_truth_Levine_BMMC_32.pdf", 
                          width = 7, height = 4.5)

# extract column order from heatmap hierarchical clustering

cols_order <- heatmap_truth$tree_col$order
cols_order

heatmap_truth$tree_col$labels  # alternatively use labels



# ========================
# cluster medians: ACCENSE
# ========================

# note ACCENSE used subsampling, so need to use data matrix with subsampled points only

cols_remove <- c("file", "label", "SNEx", "SNEy", "population")

data_medians_ACCENSE <- as.data.frame(data_ACCENSE[, -match(cols_remove, colnames(data_ACCENSE))])
dim(data_medians_ACCENSE)

colnames(data_medians_ACCENSE) <- gsub("^X\\.([A-Za-z0-9]+)\\..*$", "\\1", colnames(data_medians_ACCENSE))

table(clus_ACCENSE)
length(clus_ACCENSE)

# calculate cluster medians

medians_ACCENSE <- calculate_cluster_medians_scaled(data_medians_ACCENSE, clus_ACCENSE)

# re-order columns

medians_ACCENSE <- medians_ACCENSE[, cols_order]

# plot

pheatmap(medians_ACCENSE, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: ACCENSE", 
         filename = "../plots/cluster_medians_ACCENSE_Levine_BMMC_32.pdf", 
         width = 7, height = 7)



# =======================
# cluster medians: DensVM
# =======================

# note DensVM used subsampling, so need to use data matrix with subsampled points only

cols_remove <- c("label", "tsne_1", "tsne_2", "clust_cor_1", "clust_cor_2")

file_truth_DensVM <- "../results/Levine_BMMC_32/DensVM/cytofkit_analysis_analyzedFCS/Levine_BMMC_32_notransf.fcs"
data_truth_DensVM <- flowCore::exprs(flowCore::read.FCS(file_truth_DensVM, transformation = FALSE))

data_medians_DensVM <- data_truth_DensVM[, -match(cols_remove, colnames(data_truth_DensVM))]
dim(data_medians_DensVM)

colnames(data_medians_DensVM) <- gsub("\\(.*$", "", colnames(data_medians_DensVM))

table(clus_DensVM)
length(clus_DensVM)

# calculate cluster medians

medians_DensVM <- calculate_cluster_medians_scaled(data_medians_DensVM, clus_DensVM)

# re-order columns

medians_DensVM <- medians_DensVM[, cols_order]

# plot

pheatmap(medians_DensVM, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: DensVM", 
         filename = "../plots/cluster_medians_DensVM_Levine_BMMC_32.pdf", 
         width = 7, height = 3.5)



# ==========================
# cluster medians: flowMeans
# ==========================

dim(data_medians)

table(clus_flowMeans)
length(clus_flowMeans)

# calculate cluster medians

medians_flowMeans <- calculate_cluster_medians_scaled(data_medians, clus_flowMeans)

# re-order columns

medians_flowMeans <- medians_flowMeans[, cols_order]

# plot

pheatmap(medians_flowMeans, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: flowMeans", 
         filename = "../plots/cluster_medians_flowMeans_Levine_BMMC_32.pdf", 
         width = 7, height = 4.5)



# ========================
# cluster medians: FlowSOM
# ========================

dim(data_medians)

table(clus_FlowSOM)
length(clus_FlowSOM)

# calculate cluster medians

medians_FlowSOM <- calculate_cluster_medians_scaled(data_medians, clus_FlowSOM)

# re-order columns

medians_FlowSOM <- medians_FlowSOM[, cols_order]

# plot

pheatmap(medians_FlowSOM, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: FlowSOM", 
         filename = "../plots/cluster_medians_FlowSOM_Levine_BMMC_32.pdf", 
         width = 7, height = 12)



# =============================
# cluster medians: FlowSOM_meta
# =============================

dim(data_medians)

table(clus_FlowSOM_meta)
length(clus_FlowSOM_meta)

# calculate cluster medians

medians_FlowSOM_meta <- calculate_cluster_medians_scaled(data_medians, clus_FlowSOM_meta)

# re-order columns

medians_FlowSOM_meta <- medians_FlowSOM_meta[, cols_order]

# plot

pheatmap(medians_FlowSOM_meta, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: FlowSOM_meta", 
         filename = "../plots/cluster_medians_FlowSOM_meta_Levine_BMMC_32.pdf", 
         width = 7, height = 5)



# ============================
# cluster medians: immunoClust
# ============================

dim(data_medians)

table(clus_immunoClust)
length(clus_immunoClust)

# calculate cluster medians

medians_immunoClust <- calculate_cluster_medians_scaled(data_medians, clus_immunoClust)

# re-order columns

medians_immunoClust <- medians_immunoClust[, cols_order]

# plot

pheatmap(medians_immunoClust, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: immunoClust", 
         filename = "../plots/cluster_medians_immunoClust_Levine_BMMC_32.pdf", 
         width = 7, height = 10)



# ================================
# cluster medians: immunoClust_all
# ================================

dim(data_medians)

table(clus_immunoClust_all)
length(clus_immunoClust_all)

# calculate cluster medians

medians_immunoClust_all <- calculate_cluster_medians_scaled(data_medians, clus_immunoClust_all)

# re-order columns

medians_immunoClust_all <- medians_immunoClust_all[, cols_order]

# plot

pheatmap(medians_immunoClust_all, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: immunoClust_all", 
         filename = "../plots/cluster_medians_immunoClust_all_Levine_BMMC_32.pdf", 
         width = 7, height = 10)



# =======================
# cluster medians: kmeans
# =======================

dim(data_medians)

table(clus_kmeans)
length(clus_kmeans)

# calculate cluster medians

medians_kmeans <- calculate_cluster_medians_scaled(data_medians, clus_kmeans)

# re-order columns

medians_kmeans <- medians_kmeans[, cols_order]

# plot

pheatmap(medians_kmeans, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: kmeans", 
         filename = "../plots/cluster_medians_kmeans_Levine_BMMC_32.pdf", 
         width = 7, height = 4)



# ===========================
# cluster medians: PhenoGraph
# ===========================

dim(data_medians)

table(clus_PhenoGraph)
length(clus_PhenoGraph)

# calculate cluster medians

medians_PhenoGraph <- calculate_cluster_medians_scaled(data_medians, clus_PhenoGraph)

# re-order columns

medians_PhenoGraph <- medians_PhenoGraph[, cols_order]

# plot

pheatmap(medians_PhenoGraph, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: PhenoGraph", 
         filename = "../plots/cluster_medians_PhenoGraph_Levine_BMMC_32.pdf", 
         width = 7, height = 5)



# ===========================
# cluster medians: Rclusterpp
# ===========================

dim(data_medians)

table(clus_Rclusterpp)
length(clus_Rclusterpp)

# calculate cluster medians

medians_Rclusterpp <- calculate_cluster_medians_scaled(data_medians, clus_Rclusterpp)

# re-order columns

medians_Rclusterpp <- medians_Rclusterpp[, cols_order]

# plot

pheatmap(medians_Rclusterpp, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: Rclusterpp", 
         filename = "../plots/cluster_medians_Rclusterpp_Levine_BMMC_32.pdf", 
         width = 7, height = 5)



# ============================
# cluster medians: SamSPECTRAL
# ============================

dim(data_medians)

table(clus_SamSPECTRAL)
length(clus_SamSPECTRAL)

# calculate cluster medians

medians_SamSPECTRAL <- calculate_cluster_medians_scaled(data_medians, clus_SamSPECTRAL)

# re-order columns

medians_SamSPECTRAL <- medians_SamSPECTRAL[, cols_order]

# plot

pheatmap(medians_SamSPECTRAL, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         main = "Cluster medians: SamSPECTRAL", 
         filename = "../plots/cluster_medians_SamSPECTRAL_Levine_BMMC_32.pdf", 
         width = 7, height = 3)



# ======================
# cluster medians: SWIFT
# ======================

dim(data_medians)

table(clus_SWIFT)
length(clus_SWIFT)

# calculate cluster medians

medians_SWIFT <- calculate_cluster_medians_scaled(data_medians, clus_SWIFT)

# re-order columns

medians_SWIFT <- medians_SWIFT[, cols_order]

# plot

pheatmap(medians_SWIFT, 
         color = colorRampPalette(brewer.pal(7, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = FALSE,  # note no column clustering
         labels_row = rep("", nrow(medians_SWIFT)), 
         main = "Cluster medians: SWIFT", 
         filename = "../plots/cluster_medians_SWIFT_Levine_BMMC_32.pdf", 
         width = 7, height = 12)


