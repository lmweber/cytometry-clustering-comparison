#########################################################################################
# R script to calculate and plot median marker expression for each cluster, and compare 
# against manually gated (truth)
#
# Lukas M. Weber, November 2015
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

# data frame containing expression columns only

marker_cols_Levine <- 5:36
data_medians_Levine <- data_truth_Levine[, marker_cols_Levine]

# list of cluster labels

clus_Levine <- list(ACCENSE = clus_ACCENSE_Levine, 
                    DensVM = clus_DensVM_Levine, 
                    FLOCK = clus_FLOCK_Levine, 
                    flowMeans = clus_flowMeans_Levine, 
                    FlowSOM = clus_FlowSOM_Levine, 
                    FlowSOM_meta = clus_FlowSOM_meta_Levine, 
                    immunoClust = clus_immunoClust_Levine, 
                    immunoClust_all = clus_immunoClust_all_Levine, 
                    kmeans = clus_kmeans_Levine, 
                    PhenoGraph = clus_PhenoGraph_Levine, 
                    Rclusterpp = clus_Rclusterpp_Levine, 
                    SamSPECTRAL = clus_SamSPECTRAL_Levine, 
                    SWIFT = clus_SWIFT_Levine)

n_methods <- length(clus_Levine)



####################################################
### PLOT CLUSTER MEDIANS: MANUALLY GATED (TRUTH) ###
####################################################

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_truth_Levine <- helper_calculate_cluster_medians(data_medians_Levine, clus_truth_Levine)

rownames(medians_truth_Levine) <- paste0("manually_gated_", rownames(medians_truth_Levine))

# plot

set.seed(123)
pheatmap(medians_truth_Levine, 
         color = colorRampPalette(brewer.pal(9, "Oranges"))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Cluster medians: manually gated", 
         filename = "../plots/cluster_medians_truth_Levine_2015_marrow_32.pdf", 
         width = 8, 
         height = 4.25)



###################################################
### PLOT CLUSTER MEDIANS: EACH METHOD VS. TRUTH ###
###################################################

medians_Levine <- list()


# -------------------------------------
# calculate cluster medians for ACCENSE
# -------------------------------------

# note ACCENSE used subsampling, so need to use data matrix with subsampled points only

marker_cols_ACCENSE_Levine <- 6:37
data_medians_ACCENSE_Levine <- data_ACCENSE_Levine[, marker_cols_ACCENSE_Levine]

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_ACCENSE_Levine <- helper_calculate_cluster_medians(data_medians_ACCENSE_Levine, clus_ACCENSE_Levine)
rownames(medians_ACCENSE_Levine) <- paste0("ACCENSE_", rownames(medians_ACCENSE_Levine))

medians_Levine[[1]] <- medians_ACCENSE_Levine


# ------------------------------------
# calculate cluster medians for DensVM
# ------------------------------------

# note DensVM used subsampling, so need to use data matrix with subsampled points only

file_DensVM_Levine_sub <- "../results/DensVM/Levine_2015_marrow_32/cytofkit_analysis_analyzedFCS/Levine_2015_marrow_32_notransf.fcs"
data_medians_DensVM_Levine <- flowCore::exprs(flowCore::read.FCS(file_DensVM_Levine_sub, transformation = FALSE))

marker_cols_DensVM_Levine <- 5:36
data_medians_DensVM_Levine <- data_medians_DensVM_Levine[, marker_cols_DensVM_Levine]

# calculate cluster medians
# note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1

medians_DensVM_Levine <- helper_calculate_cluster_medians(data_medians_DensVM_Levine, clus_DensVM_Levine)
rownames(medians_DensVM_Levine) <- paste0("DensVM_", rownames(medians_DensVM_Levine))

medians_Levine[[2]] <- medians_DensVM_Levine


# -----------------------------------------------
# calculate cluster medians for all other methods
# -----------------------------------------------

for (i in 3:n_methods) {
  # calculate cluster medians
  # note values are already asinh transformed, and each dimension will be scaled to min = 0, max = 1
  
  medians_i <- helper_calculate_cluster_medians(data_medians_Levine, clus_Levine[[i]])
  rownames(medians_i) <- paste0(names(clus_Levine)[i], "_", rownames(medians_i))
  
  medians_Levine[[i]] <- medians_i
}

# method names
names(medians_Levine) <- names(clus_Levine)


# --------------
# generate plots
# --------------

plot_heights <- c(9.5, 7, 8, 8, 14, 8, 14, 14, 8, 8.5, 8, 6.5, 14)
fontsize_row <- c(10, 10, 10, 10, 7, 10, 8, 8, 10, 10, 10, 10, 2)

# plot each method together with manually gated clusters

for (i in 1:n_methods) {
  
  data_heatmap <- rbind(medians_truth_Levine, medians_Levine[[i]])
  
  annot_row <- data.frame(method = rep(c("manually_gated", names(medians_Levine)[i]), 
                                       times = c(nrow(medians_truth_Levine), nrow(medians_Levine[[i]]))))
  rownames(annot_row) <- rownames(data_heatmap)
  
  annot_colors <- c("red", "blue")
  names(annot_colors) <- c("manually_gated", names(medians_Levine)[i])
  annot_colors <- list(method = annot_colors)
  
  set.seed(123)
  pheatmap(data_heatmap, 
           color = colorRampPalette(brewer.pal(9, "YlGn"))(100), 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           annotation_row = annot_row, 
           annotation_colors = annot_colors, 
           fontsize_row = fontsize_row[i], 
           main = paste0("Cluster medians: manually gated and ", names(medians_Levine)[i]), 
           filename = paste0("../plots/cluster_medians_", names(medians_Levine)[i], "_Levine_2015_marrow_32.pdf"), 
           width = 10, 
           height = plot_heights[i])
}


