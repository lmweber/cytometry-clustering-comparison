#########################################################################################
# R script to calculate and plot median marker expression for each cluster, and compare 
# against manually gated (truth)
#
# Lukas M. Weber, December 2015
#########################################################################################


library(pheatmap)
library(RColorBrewer)

# helper function
source("helper_calculate_cluster_medians.R")

# load results from previous steps
source("load_results_ACCENSE.R")
source("load_results_DensVM.R")
source("load_results_FLOCK.R")
source("load_results_PhenoGraph.R")
source("load_results_SWIFT.R")
source("load_results_truth.R")
source("load_results_all_other_methods.R")




######################
### DATA FOR PLOTS ###
######################

# data frames containing expression columns only

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13

data_medians_Levine_32 <- data_truth_Levine_32[, marker_cols_Levine_32]
data_medians_Levine_13 <- data_truth_Levine_13[, marker_cols_Levine_13]


# lists of cluster labels

clus_Levine_32 <- list(ACCENSE = clus_ACCENSE_Levine_32, 
                       DensVM = clus_DensVM_Levine_32, 
                       FLOCK = clus_FLOCK_Levine_32, 
                       flowMeans = clus_flowMeans_Levine_32, 
                       FlowSOM = clus_FlowSOM_Levine_32, 
                       FlowSOM_meta = clus_FlowSOM_meta_Levine_32, 
                       immunoClust = clus_immunoClust_Levine_32, 
                       immunoClust_all = clus_immunoClust_all_Levine_32, 
                       kmeans = clus_kmeans_Levine_32, 
                       PhenoGraph = clus_PhenoGraph_Levine_32, 
                       Rclusterpp = clus_Rclusterpp_Levine_32, 
                       SamSPECTRAL = clus_SamSPECTRAL_Levine_32, 
                       SWIFT = clus_SWIFT_Levine_32)

clus_Levine_13 <- list(ACCENSE = clus_ACCENSE_Levine_13, 
                       DensVM = clus_DensVM_Levine_13, 
                       FLOCK = clus_FLOCK_Levine_13, 
                       flowMeans = clus_flowMeans_Levine_13, 
                       FlowSOM = clus_FlowSOM_Levine_13, 
                       FlowSOM_meta = clus_FlowSOM_meta_Levine_13, 
                       immunoClust = clus_immunoClust_Levine_13, 
                       immunoClust_all = clus_immunoClust_all_Levine_13, 
                       kmeans = clus_kmeans_Levine_13, 
                       PhenoGraph = clus_PhenoGraph_Levine_13, 
                       Rclusterpp = clus_Rclusterpp_Levine_13, 
                       SamSPECTRAL = clus_SamSPECTRAL_Levine_13, 
                       SWIFT = clus_SWIFT_Levine_13)


n_methods_Levine_32 <- length(clus_Levine_32)
n_methods_Levine_13 <- length(clus_Levine_13)




####################################################
### PLOT CLUSTER MEDIANS: MANUALLY GATED (TRUTH) ###
####################################################

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_truth_Levine_32 <- helper_calculate_cluster_medians(data_medians_Levine_32, clus_truth_Levine_32)
medians_truth_Levine_13 <- helper_calculate_cluster_medians(data_medians_Levine_13, clus_truth_Levine_13)

rownames(medians_truth_Levine_32) <- paste0("manually_gated_", rownames(medians_truth_Levine_32))
rownames(medians_truth_Levine_13) <- paste0("manually_gated_", rownames(medians_truth_Levine_13))


# plot heatmaps

set.seed(123)
pheatmap(medians_truth_Levine_32, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         fontsize = 9, 
         filename = "../plots/Levine_2015_marrow_32/cluster_medians_truth_Levine2015marrow32.pdf", 
         width = 8, 
         height = 4)

set.seed(123)
pheatmap(medians_truth_Levine_13, 
         color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         fontsize = 9, 
         filename = "../plots/Levine_2015_marrow_13/cluster_medians_truth_Levine2015marrow13.pdf", 
         width = 5, 
         height = 6)




###########################################################
### PLOT CLUSTER MEDIANS: EACH METHOD COMPARED TO TRUTH ###
###########################################################

medians_Levine_32 <- list()
medians_Levine_13 <- list()


### calculate cluster medians for ACCENSE ###

# note ACCENSE used subsampling, so need to use data matrix with subsampled points only

marker_cols_ACCENSE_Levine_32 <- 6:37
marker_cols_ACCENSE_Levine_13 <- 2:14

data_medians_ACCENSE_Levine_32 <- data_ACCENSE_Levine_32[, marker_cols_ACCENSE_Levine_32]
data_medians_ACCENSE_Levine_13 <- data_ACCENSE_Levine_13[, marker_cols_ACCENSE_Levine_13]

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_ACCENSE_Levine_32 <- helper_calculate_cluster_medians(data_medians_ACCENSE_Levine_32, clus_ACCENSE_Levine_32)
medians_ACCENSE_Levine_13 <- helper_calculate_cluster_medians(data_medians_ACCENSE_Levine_13, clus_ACCENSE_Levine_13)

rownames(medians_ACCENSE_Levine_32) <- paste0("ACCENSE_", rownames(medians_ACCENSE_Levine_32))
rownames(medians_ACCENSE_Levine_13) <- paste0("ACCENSE_", rownames(medians_ACCENSE_Levine_13))

medians_Levine_32[[1]] <- medians_ACCENSE_Levine_32
medians_Levine_13[[1]] <- medians_ACCENSE_Levine_13


### calculate cluster medians for DensVM ###

# note DensVM used subsampling, so need to use data matrix with subsampled points only

file_DensVM_Levine_32_sub <- "../results/DensVM/Levine_2015_marrow_32/cytofkit_analysis_analyzedFCS/Levine_2015_marrow_32_notransform.fcs"
file_DensVM_Levine_13_sub <- "../results/DensVM/Levine_2015_marrow_13/cytofkit_analysis_analyzedFCS/Levine_2015_marrow_13_notransform.fcs"

data_medians_DensVM_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_DensVM_Levine_32_sub, transformation = FALSE))
data_medians_DensVM_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_DensVM_Levine_13_sub, transformation = FALSE))

marker_cols_DensVM_Levine_32 <- 5:36
marker_cols_DensVM_Levine_13 <- 1:13

data_medians_DensVM_Levine_32 <- data_medians_DensVM_Levine_32[, marker_cols_DensVM_Levine_32]
data_medians_DensVM_Levine_13 <- data_medians_DensVM_Levine_13[, marker_cols_DensVM_Levine_13]

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_DensVM_Levine_32 <- helper_calculate_cluster_medians(data_medians_DensVM_Levine_32, clus_DensVM_Levine_32)
medians_DensVM_Levine_13 <- helper_calculate_cluster_medians(data_medians_DensVM_Levine_13, clus_DensVM_Levine_13)

rownames(medians_DensVM_Levine_32) <- paste0("DensVM_", rownames(medians_DensVM_Levine_32))
rownames(medians_DensVM_Levine_13) <- paste0("DensVM_", rownames(medians_DensVM_Levine_13))

medians_Levine_32[[2]] <- medians_DensVM_Levine_32
medians_Levine_13[[2]] <- medians_DensVM_Levine_13


### calculate cluster medians for all other methods ###

for (i in 3:n_methods_Levine_32) {
  # calculate cluster medians
  # note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1
  medians_i <- helper_calculate_cluster_medians(data_medians_Levine_32, clus_Levine_32[[i]])
  rownames(medians_i) <- paste0(names(clus_Levine_32)[i], "_", rownames(medians_i))
  medians_Levine_32[[i]] <- medians_i
}

for (i in 3:n_methods_Levine_13) {
  # calculate cluster medians
  # note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1
  medians_i <- helper_calculate_cluster_medians(data_medians_Levine_13, clus_Levine_13[[i]])
  rownames(medians_i) <- paste0(names(clus_Levine_13)[i], "_", rownames(medians_i))
  medians_Levine_13[[i]] <- medians_i
}

# method names

names(medians_Levine_32) <- names(clus_Levine_32)
names(medians_Levine_13) <- names(clus_Levine_13)


### generate plots ###

# plot each method together with manually gated clusters

plot_heights_Levine_32 <- c(8.5, 6.5, 7.5, 7, 14, 7.5, 14, 14, 7.5, 8, 7.5, 6, 14)
plot_heights_Levine_13 <- c(10.5, 8.5, 9.5, 9, 14, 12, 14, 14, 12, 9, 12, 8, 14)

fontsize_row_Levine_32 <- c(9, 9, 9, 9, 7, 9, 8, 8, 9, 9, 9, 9, 2)
fontsize_row_Levine_13 <- c(9, 9, 9, 9, 7, 9, 7, 7, 9, 9, 9, 9, 4)


for (i in 1:n_methods_Levine_32) {
  
  data_heatmap <- rbind(medians_truth_Levine_32, medians_Levine_32[[i]])
  
  annot_row <- data.frame(method = rep(c("manually_gated", names(medians_Levine_32)[i]), 
                                       times = c(nrow(medians_truth_Levine_32), nrow(medians_Levine_32[[i]]))))
  rownames(annot_row) <- rownames(data_heatmap)
  
  annot_colors <- c("red", "blue")
  names(annot_colors) <- c("manually_gated", names(medians_Levine_32)[i])
  annot_colors <- list(method = annot_colors)
  
  set.seed(123)
  pheatmap(data_heatmap, 
           color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           annotation_row = annot_row, 
           annotation_colors = annot_colors, 
           fontsize = 9, 
           fontsize_row = fontsize_row_Levine_32[i], 
           filename = paste0("../plots/Levine_2015_marrow_32/cluster_medians_", names(medians_Levine_32)[i], "_Levine2015marrow32.pdf"), 
           width = 10, 
           height = plot_heights_Levine_32[i])
}


for (i in 1:n_methods_Levine_13) {
  
  data_heatmap <- rbind(medians_truth_Levine_13, medians_Levine_13[[i]])
  
  annot_row <- data.frame(method = rep(c("manually_gated", names(medians_Levine_13)[i]), 
                                       times = c(nrow(medians_truth_Levine_13), nrow(medians_Levine_13[[i]]))))
  rownames(annot_row) <- rownames(data_heatmap)
  
  annot_colors <- c("red", "blue")
  names(annot_colors) <- c("manually_gated", names(medians_Levine_13)[i])
  annot_colors <- list(method = annot_colors)
  
  set.seed(123)
  pheatmap(data_heatmap, 
           color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           annotation_row = annot_row, 
           annotation_colors = annot_colors, 
           fontsize = 9, 
           fontsize_row = fontsize_row_Levine_13[i], 
           filename = paste0("../plots/Levine_2015_marrow_13/cluster_medians_", names(medians_Levine_13)[i], "_Levine2015marrow13.pdf"), 
           width = 6.75, 
           height = plot_heights_Levine_13[i])
}


