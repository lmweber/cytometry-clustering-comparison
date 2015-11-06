#########################################################################################
# R script to identify clusters: calculate median marker expression for each cluster, 
# then use heatmaps with hierarchical clustering to compare with manual gating results
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

source("calculate_cluster_medians_scaled.R")

# load results from previous steps

source("heatmaps_tables_Levine_BMMC_32.R")



#########################################################################################
### Calculate cluster medians and compare each method with manually gated populations ###
#########################################################################################


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

rownames(medians_truth) <- paste0("manually_gated_", rownames(medians_truth))

# plot

pheatmap(medians_truth, 
         color = colorRampPalette(brewer.pal(9, "Oranges"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         main = "Cluster medians: manually gated", 
         filename = "../plots/cluster_medians_truth_Levine_BMMC_32.pdf", 
         width = 8, height = 4.25)



# ========================
# cluster medians: ACCENSE
# ========================

# note ACCENSE used subsampling, so need to use data matrix with subsampled points only

cols_remove <- c("file", "label", "SNEx", "SNEy", "population")

data_medians_ACCENSE <- as.data.frame(data_ACCENSE[, -match(cols_remove, colnames(data_ACCENSE))])
dim(data_medians_ACCENSE)

colnames(data_medians_ACCENSE) <- gsub("(^X\\.)|(\\.[A-Za-z0-9]+\\.Di$)", "", colnames(data_medians_ACCENSE))

table(clus_ACCENSE)
length(clus_ACCENSE)

# calculate cluster medians

medians_ACCENSE <- calculate_cluster_medians_scaled(data_medians_ACCENSE, clus_ACCENSE)

rownames(medians_ACCENSE) <- paste0("ACCENSE_", rownames(medians_ACCENSE))

# plot

data_heatmap <- rbind(medians_truth, medians_ACCENSE)
annot_row <- data.frame(method = rep(c("manually_gated", "ACCENSE"), 
                                     times = c(nrow(medians_truth), nrow(medians_ACCENSE))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", ACCENSE = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and ACCENSE", 
         filename = "../plots/cluster_medians_ACCENSE_Levine_BMMC_32.pdf", 
         width = 10, height = 10)



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

rownames(medians_DensVM) <- paste0("DensVM_", rownames(medians_DensVM))

# plot

data_heatmap <- rbind(medians_truth, medians_DensVM)
annot_row <- data.frame(method = rep(c("manually_gated", "DensVM"), 
                                     times = c(nrow(medians_truth), nrow(medians_DensVM))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", DensVM = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and DensVM", 
         filename = "../plots/cluster_medians_DensVM_Levine_BMMC_32.pdf", 
         width = 10, height = 7)



# ==========================
# cluster medians: flowMeans
# ==========================

dim(data_medians)

table(clus_flowMeans)
length(clus_flowMeans)

# calculate cluster medians

medians_flowMeans <- calculate_cluster_medians_scaled(data_medians, clus_flowMeans)

rownames(medians_flowMeans) <- paste0("flowMeans_", rownames(medians_flowMeans))

# plot

data_heatmap <- rbind(medians_truth, medians_flowMeans)
annot_row <- data.frame(method = rep(c("manually_gated", "flowMeans"), 
                                     times = c(nrow(medians_truth), nrow(medians_flowMeans))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", flowMeans = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and flowMeans", 
         filename = "../plots/cluster_medians_flowMeans_Levine_BMMC_32.pdf", 
         width = 10, height = 8)



# ========================
# cluster medians: FlowSOM
# ========================

dim(data_medians)

table(clus_FlowSOM)
length(clus_FlowSOM)

# calculate cluster medians

medians_FlowSOM <- calculate_cluster_medians_scaled(data_medians, clus_FlowSOM)

rownames(medians_FlowSOM) <- paste0("FlowSOM_", rownames(medians_FlowSOM))

# plot

data_heatmap <- rbind(medians_truth, medians_FlowSOM)
annot_row <- data.frame(method = rep(c("manually_gated", "FlowSOM"), 
                                     times = c(nrow(medians_truth), nrow(medians_FlowSOM))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", FlowSOM = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         fontsize_row = 7, 
         main = "Cluster medians: manually gated and FlowSOM", 
         filename = "../plots/cluster_medians_FlowSOM_Levine_BMMC_32.pdf", 
         width = 10, height = 14)



# =============================
# cluster medians: FlowSOM_meta
# =============================

dim(data_medians)

table(clus_FlowSOM_meta)
length(clus_FlowSOM_meta)

# calculate cluster medians

medians_FlowSOM_meta <- calculate_cluster_medians_scaled(data_medians, clus_FlowSOM_meta)

rownames(medians_FlowSOM_meta) <- paste0("FlowSOM_meta_", rownames(medians_FlowSOM_meta))

# plot

data_heatmap <- rbind(medians_truth, medians_FlowSOM_meta)
annot_row <- data.frame(method = rep(c("manually_gated", "FlowSOM_meta"), 
                                     times = c(nrow(medians_truth), nrow(medians_FlowSOM_meta))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", FlowSOM_meta = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and FlowSOM_meta", 
         filename = "../plots/cluster_medians_FlowSOM_meta_Levine_BMMC_32.pdf", 
         width = 10, height = 8)



# ============================
# cluster medians: immunoClust
# ============================

dim(data_medians)

table(clus_immunoClust)
length(clus_immunoClust)

# calculate cluster medians

medians_immunoClust <- calculate_cluster_medians_scaled(data_medians, clus_immunoClust)

rownames(medians_immunoClust) <- paste0("immunoClust_", rownames(medians_immunoClust))

# plot

data_heatmap <- rbind(medians_truth, medians_immunoClust)
annot_row <- data.frame(method = rep(c("manually_gated", "immunoClust"), 
                                     times = c(nrow(medians_truth), nrow(medians_immunoClust))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", immunoClust = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         fontsize_row = 8, 
         main = "Cluster medians: manually gated and immunoClust", 
         filename = "../plots/cluster_medians_immunoClust_Levine_BMMC_32.pdf", 
         width = 10, height = 14)



# ================================
# cluster medians: immunoClust_all
# ================================

dim(data_medians)

table(clus_immunoClust_all)
length(clus_immunoClust_all)

# calculate cluster medians

medians_immunoClust_all <- calculate_cluster_medians_scaled(data_medians, clus_immunoClust_all)

rownames(medians_immunoClust_all) <- paste0("immunoClust_all_", rownames(medians_immunoClust_all))

# plot

data_heatmap <- rbind(medians_truth, medians_immunoClust_all)
annot_row <- data.frame(method = rep(c("manually_gated", "immunoClust_all"), 
                                     times = c(nrow(medians_truth), nrow(medians_immunoClust_all))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", immunoClust_all = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         fontsize_row = 8, 
         main = "Cluster medians: manually gated and immunoClust_all", 
         filename = "../plots/cluster_medians_immunoClust_all_Levine_BMMC_32.pdf", 
         width = 10, height = 14)



# =======================
# cluster medians: kmeans
# =======================

dim(data_medians)

table(clus_kmeans)
length(clus_kmeans)

# calculate cluster medians

medians_kmeans <- calculate_cluster_medians_scaled(data_medians, clus_kmeans)

rownames(medians_kmeans) <- paste0("kmeans_", rownames(medians_kmeans))

# plot

data_heatmap <- rbind(medians_truth, medians_kmeans)
annot_row <- data.frame(method = rep(c("manually_gated", "kmeans"), 
                                     times = c(nrow(medians_truth), nrow(medians_kmeans))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", kmeans = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and kmeans", 
         filename = "../plots/cluster_medians_kmeans_Levine_BMMC_32.pdf", 
         width = 10, height = 8)



# ===========================
# cluster medians: PhenoGraph
# ===========================

dim(data_medians)

table(clus_PhenoGraph)
length(clus_PhenoGraph)

# calculate cluster medians

medians_PhenoGraph <- calculate_cluster_medians_scaled(data_medians, clus_PhenoGraph)

rownames(medians_PhenoGraph) <- paste0("PhenoGraph_", rownames(medians_PhenoGraph))

# plot

data_heatmap <- rbind(medians_truth, medians_PhenoGraph)
annot_row <- data.frame(method = rep(c("manually_gated", "PhenoGraph"), 
                                     times = c(nrow(medians_truth), nrow(medians_PhenoGraph))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", PhenoGraph = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and PhenoGraph", 
         filename = "../plots/cluster_medians_PhenoGraph_Levine_BMMC_32.pdf", 
         width = 10, height = 8.5)



# ===========================
# cluster medians: Rclusterpp
# ===========================

dim(data_medians)

table(clus_Rclusterpp)
length(clus_Rclusterpp)

# calculate cluster medians

medians_Rclusterpp <- calculate_cluster_medians_scaled(data_medians, clus_Rclusterpp)

rownames(medians_Rclusterpp) <- paste0("Rclusterpp_", rownames(medians_Rclusterpp))

# plot

data_heatmap <- rbind(medians_truth, medians_Rclusterpp)
annot_row <- data.frame(method = rep(c("manually_gated", "Rclusterpp"), 
                                     times = c(nrow(medians_truth), nrow(medians_Rclusterpp))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", Rclusterpp = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and Rclusterpp", 
         filename = "../plots/cluster_medians_Rclusterpp_Levine_BMMC_32.pdf", 
         width = 10, height = 8)



# ============================
# cluster medians: SamSPECTRAL
# ============================

dim(data_medians)

table(clus_SamSPECTRAL)
length(clus_SamSPECTRAL)

# calculate cluster medians

medians_SamSPECTRAL <- calculate_cluster_medians_scaled(data_medians, clus_SamSPECTRAL)

rownames(medians_SamSPECTRAL) <- paste0("SamSPECTRAL_", rownames(medians_SamSPECTRAL))

# plot

data_heatmap <- rbind(medians_truth, medians_SamSPECTRAL)
annot_row <- data.frame(method = rep(c("manually_gated", "SamSPECTRAL"), 
                                     times = c(nrow(medians_truth), nrow(medians_SamSPECTRAL))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", SamSPECTRAL = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         main = "Cluster medians: manually gated and SamSPECTRAL", 
         filename = "../plots/cluster_medians_SamSPECTRAL_Levine_BMMC_32.pdf", 
         width = 10, height = 6.5)



# ======================
# cluster medians: SWIFT
# ======================

dim(data_medians)

table(clus_SWIFT)
length(clus_SWIFT)

# calculate cluster medians

medians_SWIFT <- calculate_cluster_medians_scaled(data_medians, clus_SWIFT)

rownames(medians_SWIFT) <- paste0("SWIFT_", rownames(medians_SWIFT))

# plot

data_heatmap <- rbind(medians_truth, medians_SWIFT)
annot_row <- data.frame(method = rep(c("manually_gated", "SWIFT"), 
                                     times = c(nrow(medians_truth), nrow(medians_SWIFT))))
rownames(annot_row) <- rownames(data_heatmap)
annot_colors <- list(method = c(manually_gated = "red", SWIFT = "blue"))

pheatmap(data_heatmap, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_row = annot_row, annotation_colors = annot_colors, 
         fontsize_row = 2, 
         main = "Cluster medians: manually gated and SWIFT", 
         filename = "../plots/cluster_medians_SWIFT_Levine_BMMC_32.pdf", 
         width = 10, height = 14)


