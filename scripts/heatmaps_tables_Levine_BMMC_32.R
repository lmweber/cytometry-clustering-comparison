#########################################################################################
# R script to generate heatmaps and tables of results comparing clustering methods
#
# data set "Levine_BMMC_32"
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)
library(pheatmap)
library(RColorBrewer)

source("match_clusters_and_evaluate.R")

source("load_results_Levine_BMMC_32.R")



################
### HEATMAPS ###
################


# ======================================================================
# Plot heatmaps of precision, recall, and F1 score for each true cluster
# ======================================================================

res_pr <- cbind(pr_ACCENSE, pr_DensVM, pr_FlowSOM, pr_FlowSOM_meta, pr_immunoClust, pr_immunoClust_all, 
                pr_kmeans, pr_PhenoGraph, pr_Rclusterpp, pr_SWIFT)
res_re <- cbind(re_ACCENSE, re_DensVM, re_FlowSOM, re_FlowSOM_meta, re_immunoClust, re_immunoClust_all, 
                re_kmeans, re_PhenoGraph, re_Rclusterpp, re_SWIFT)
res_F1 <- cbind(F1_ACCENSE, F1_DensVM, F1_FlowSOM, F1_FlowSOM_meta, F1_immunoClust, F1_immunoClust_all, 
                F1_kmeans, F1_PhenoGraph, F1_Rclusterpp, F1_SWIFT)

colnames(res_pr) <- gsub("^pr_", "", colnames(res_pr))
colnames(res_re) <- gsub("^re_", "", colnames(res_re))
colnames(res_F1) <- gsub("^F1_", "", colnames(res_F1))

png("../plots/heatmap_precision_Levine_BMMC_32.png", width = 1000, height = 1200)
pheatmap(res_pr, 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
         breaks = seq(0, 1, length.out = 100), 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, fontsize_number = 13, cex = 1.5, 
         main = "precision: Levine_BMMC_32")
dev.off()

png("../plots/heatmap_recall_Levine_BMMC_32.png", width = 1000, height = 1200)
pheatmap(res_re, 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
         breaks = seq(0, 1, length.out = 100), 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, fontsize_number = 13, cex = 1.5, 
         main = "recall: Levine_BMMC_32")
dev.off()

png("../plots/heatmap_F1_Levine_BMMC_32.png", width = 1000, height = 1200)
pheatmap(res_F1, 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
         breaks = seq(0, 1, length.out = 100), 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, fontsize_number = 13, cex = 1.5, 
         main = "F1 score: Levine_BMMC_32")
dev.off()




#########################
### TABLES OF RESULTS ###
#########################


# ========================================
# Number of cells for each matched cluster
# ========================================

n_cells_matched <- cbind(ACCENSE = n_cells_ACCENSE, 
                         DensVM = n_cells_DensVM, 
                         FlowSOM = n_cells_FlowSOM, 
                         FlowSOM_meta = n_cells_FlowSOM_meta, 
                         immunoClust = n_cells_immunoClust, 
                         immunoClust_all = n_cells_immunoClust_all, 
                         kmeans = n_cells_kmeans, 
                         PhenoGraph = n_cells_PhenoGraph, 
                         Rclusterpp = n_cells_Rclusterpp, 
                         SWIFT = n_cells_SWIFT)

n_cells_matched_with_truth <- cbind(manually_gated = n_cells_truth, 
                                    n_cells_matched)

n_cells_matched_with_truth



# ======================
# Matched cluster labels
# ======================

labels_matched <- cbind(manually_gated = 1:n_clus_truth, 
                        ACCENSE = labels_matched_ACCENSE, 
                        DensVM = labels_matched_DensVM, 
                        FlowSOM = labels_matched_FlowSOM, 
                        FlowSOM_meta = labels_matched_FlowSOM_meta, 
                        immunoClust = labels_matched_immunoClust, 
                        immunoClust_all = labels_matched_immunoClust_all, 
                        kmeans = labels_matched_kmeans, 
                        PhenoGraph = labels_matched_PhenoGraph, 
                        Rclusterpp = labels_matched_Rclusterpp, 
                        SWIFT = labels_matched_SWIFT)
labels_matched

## note double matchings are allowed
## this is justified since just picking best match by F1 score
## and for DensVM, number of detected clusters is less than number of true clusters



# =========================================
# Number of cells for each detected cluster
# =========================================

n_cells_detected <- list(manually_gated = tbl_truth, 
                         ACCENSE = tbl_ACCENSE, 
                         DensVM = tbl_DensVM, 
                         FlowSOM = tbl_FlowSOM, 
                         FlowSOM_meta = tbl_FlowSOM_meta, 
                         immunoClust = tbl_immunoClust, 
                         immunoClust_all = tbl_immunoClust_all, 
                         kmeans = tbl_kmeans, 
                         PhenoGraph = tbl_PhenoGraph, 
                         Rclusterpp = tbl_Rclusterpp, 
                         SWIFT = tbl_SWIFT)
n_cells_detected



# ==========================================
# Number of clusters detected by each method
# ==========================================

n_clus_detected <- rbind(manually_gated = n_clus_truth, 
                         ACCENSE = n_clus_ACCENSE, 
                         DensVM = n_clus_DensVM, 
                         FlowSOM = n_clus_FlowSOM, 
                         FlowSOM_meta = n_clus_FlowSOM_meta, 
                         immunoClust = n_clus_immunoClust, 
                         immunoClust_all = n_clus_immunoClust_all, 
                         kmeans = n_clus_kmeans, 
                         PhenoGraph = n_clus_PhenoGraph, 
                         Rclusterpp = n_clus_Rclusterpp, 
                         SWIFT = n_clus_SWIFT)
n_clus_detected

# alternatively

sapply(n_cells_detected, length)



# ===============================================
# Mean F1-score (averaged over all true clusters)
# ===============================================

n_cells_truth_matrix <- matrix(rep(n_cells_truth, 10), nrow = 14)
n_cells_truth_matrix

colSums(n_cells_truth_matrix)

mean_F1 <- colSums(res_F1 * n_cells_truth_matrix) / colSums(n_cells_truth_matrix)
mean_F1


## alternatively: weight by number of matched cells

#total_cells_matched <- colSums(n_cells_matched)
#total_cells_matched

#mean_F1 <- colSums(res_F1 * n_cells_matched) / total_cells_matched
#mean_F1


png("../plots/heatmap_mean_F1_Levine_BMMC_32.png", width = 1000, height = 300)
pheatmap(t(mean_F1), 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
         breaks = seq(0, 1, length.out = 100), 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, fontsize_number = 13, cex = 1.5, 
         main = "Mean F1 score: Levine_BMMC_32")
dev.off()



# ==================
# Cluster properties
# ==================

# pull out a table of average (mean? median?) values along each marker, for each detected
# cluster as well as truth clusters
# - how similar to the true clusters are the detected clusters? do the results make
# sense, i.e. do the clusters correspond to known cell types?



