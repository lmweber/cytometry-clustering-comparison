#########################################################################################
# R script to generate main plots of results
#
# Lukas M. Weber, December 2015
#########################################################################################


library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

# helper function
source("helper_collapse_data_frame.R")

# load results from previous steps
source("load_results_ACCENSE.R")
source("load_results_DensVM.R")
source("load_results_FLOCK.R")
source("load_results_PhenoGraph.R")
source("load_results_SWIFT.R")
source("load_results_truth.R")
source("load_results_all_other_methods.R")




###########################
### PREPARE DATA FRAMES ###
###########################

# combine results into lists

res_Levine_32 <- list(ACCENSE = res_ACCENSE_Levine_32, 
                      DensVM = res_DensVM_Levine_32, 
                      FLOCK = res_FLOCK_Levine_32, 
                      flowMeans = res_flowMeans_Levine_32, 
                      FlowSOM = res_FlowSOM_Levine_32, 
                      FlowSOM_meta = res_FlowSOM_meta_Levine_32, 
                      immunoClust = res_immunoClust_Levine_32, 
                      immunoClust_all = res_immunoClust_all_Levine_32, 
                      kmeans = res_kmeans_Levine_32, 
                      PhenoGraph = res_PhenoGraph_Levine_32, 
                      Rclusterpp = res_Rclusterpp_Levine_32, 
                      SamSPECTRAL = res_SamSPECTRAL_Levine_32, 
                      SWIFT = res_SWIFT_Levine_32)

res_Levine_13 <- list(ACCENSE = res_ACCENSE_Levine_13, 
                      DensVM = res_DensVM_Levine_13, 
                      FLOCK = res_FLOCK_Levine_13, 
                      flowMeans = res_flowMeans_Levine_13, 
                      FlowSOM = res_FlowSOM_Levine_13, 
                      FlowSOM_meta = res_FlowSOM_meta_Levine_13, 
                      immunoClust = res_immunoClust_Levine_13, 
                      immunoClust_all = res_immunoClust_all_Levine_13, 
                      kmeans = res_kmeans_Levine_13, 
                      PhenoGraph = res_PhenoGraph_Levine_13, 
                      Rclusterpp = res_Rclusterpp_Levine_13, 
                      SamSPECTRAL = res_SamSPECTRAL_Levine_13, 
                      SWIFT = res_SWIFT_Levine_13)

res_Mosmann <- list(ACCENSE = res_ACCENSE_Mosmann, 
                    DensVM = res_DensVM_Mosmann, 
                    FLOCK = res_FLOCK_Mosmann, 
                    flowMeans = res_flowMeans_Mosmann, 
                    FlowSOM = res_FlowSOM_Mosmann, 
                    FlowSOM_meta = res_FlowSOM_meta_Mosmann, 
                    immunoClust = res_immunoClust_Mosmann, 
                    immunoClust_all = res_immunoClust_all_Mosmann, 
                    kmeans = res_kmeans_Mosmann, 
                    PhenoGraph = res_PhenoGraph_Mosmann, 
                    #Rclusterpp = res_Rclusterpp_Mosmann,  # Rclusterpp did not complete for this data set
                    SamSPECTRAL = res_SamSPECTRAL_Mosmann, 
                    SWIFT = res_SWIFT_Mosmann)


# collapse into data frames (using helper functions to pad missing values with zeros or NAs)

precision_df_Levine_32 <- sapply(res_Levine_32, function(res) res$pr)
precision_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$pr))
precision_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$pr))

recall_df_Levine_32 <- sapply(res_Levine_32, function(res) res$re)
recall_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$re))
recall_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$re))

F1_df_Levine_32 <- sapply(res_Levine_32, function(res) res$F1)
F1_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$F1))
F1_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$F1))

labels_matched_df_Levine_32 <- sapply(res_Levine_32, function(res) res$labels_matched)
labels_matched_df_Levine_13 <- collapse_data_frame_NAs(sapply(res_Levine_13, function(res) res$labels_matched))  # NA = no match
labels_matched_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$labels_matched))

n_cells_df_Levine_32 <- sapply(res_Levine_32, function(res) res$n_cells)
n_cells_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$n_cells))
n_cells_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$n_cells))




#####################
### MEAN F1 SCORE ###
#####################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# number of methods

n_methods_Levine_32 <- ncol(n_cells_df_Levine_32)
n_methods_Levine_13 <- ncol(n_cells_df_Levine_13)

# number of cells per true clusters

n_cells_truth_Levine_32 <- matrix(rep(tbl_truth_Levine_32, n_methods_Levine_32), ncol = n_methods_Levine_32)
n_cells_truth_Levine_13 <- matrix(rep(tbl_truth_Levine_13, n_methods_Levine_13), ncol = n_methods_Levine_13)

n_cells_truth_Levine_32
n_cells_truth_Levine_13

# mean F1 score across all true clusters, weighted by number of cells in true cluster

mean_F1_Levine_32 <- colSums(F1_df_Levine_32 * n_cells_truth_Levine_32) / colSums(n_cells_truth_Levine_32)
mean_F1_Levine_13 <- colSums(F1_df_Levine_13 * n_cells_truth_Levine_13) / colSums(n_cells_truth_Levine_13)

mean_F1_Levine_32
mean_F1_Levine_13

# arrange in descending order

ord_Levine_32 <- rev(order(mean_F1_Levine_32))
ord_Levine_13 <- rev(order(mean_F1_Levine_13))

mean_F1_Levine_32_ord <- mean_F1_Levine_32[ord_Levine_32]
mean_F1_Levine_13_ord <- mean_F1_Levine_13[ord_Levine_13]

mean_F1_Levine_32_ord
mean_F1_Levine_13_ord


# tidy data format (for ggplot)

mean_F1_Levine_32_tidy <- data.frame(value = mean_F1_Levine_32_ord)
mean_F1_Levine_32_tidy["method"] <- factor(rownames(mean_F1_Levine_32_tidy), levels = rownames(mean_F1_Levine_32_tidy))
mean_F1_Levine_32_tidy

mean_F1_Levine_13_tidy <- data.frame(value = mean_F1_Levine_13_ord)
mean_F1_Levine_13_tidy["method"] <- factor(rownames(mean_F1_Levine_13_tidy), levels = rownames(mean_F1_Levine_13_tidy))
mean_F1_Levine_13_tidy


# bar plots

ggplot(mean_F1_Levine_32_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.2f", round(value, 2)), y = value + 0.03), size = 3.5) + 
  ylim(0, 1) + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/Levine_2015_marrow_32/results_mean_F1_score_Levine2015marrow32.pdf", height = 6, width = 5.5)


ggplot(mean_F1_Levine_13_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.2f", round(value, 2)), y = value + 0.03), size = 3.5) + 
  ylim(0, 1) + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/Levine_2015_marrow_13/results_mean_F1_score_Levine2015marrow13.pdf", height = 6, width = 5.5)




########################
### POPULATION SIZES ###
########################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# tidy data format (for ggplot)

n_cells_truth_Levine_32_tidy <- as.data.frame(tbl_truth_Levine_32)
n_cells_truth_Levine_13_tidy <- as.data.frame(tbl_truth_Levine_13)

colnames(n_cells_truth_Levine_32_tidy) <- c("population", "value")
colnames(n_cells_truth_Levine_13_tidy) <- c("population", "value")


# plot number of cells in each true (manually gated) population

ggplot(n_cells_truth_Levine_32_tidy, aes(x = population, y = value)) + 
  geom_bar(stat = "identity", fill = "darkgray") + 
  geom_text(aes(label = value, y = value + 1300, angle = 90), size = 3.5) + 
  ylim(0, 28000) + 
  xlab("manually gated population") + 
  ylab("number of cells") + 
  ggtitle("Manually gated populations: Levine_2015_marrow_32") + 
  theme_bw()

ggsave("../plots/Levine_2015_marrow_32/results_number_of_cells_Levine2015marrow32.pdf", height = 6, width = 6)


ggplot(n_cells_truth_Levine_13_tidy, aes(x = population, y = value)) + 
  geom_bar(stat = "identity", fill = "darkgray") + 
  geom_text(aes(label = value, y = value + 1300, angle = 90), size = 3.5) + 
  ylim(0, 28000) + 
  xlab("manually gated population") + 
  ylab("number of cells") + 
  ggtitle("Manually gated populations: Levine_2015_marrow_13") + 
  theme_bw()

ggsave("../plots/Levine_2015_marrow_13/results_number_of_cells_Levine2015marrow13.pdf", height = 6, width = 6)




###############################################################################
### LINE PLOTS: F1 SCORE, PRECISION, AND RECALL SORTED BY TRUE CLUSTER SIZE ###
###############################################################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# order by true cluster sizes

ord_size_Levine_32 <- order(tbl_truth_Levine_32)
ord_size_Levine_13 <- order(tbl_truth_Levine_13)


# tidy data format (for ggplot)

# F1 score

F1_df_Levine_32_ord_size <- F1_df_Levine_32[ord_size_Levine_32, ]
F1_df_Levine_32_ord_size_tidy <- data.frame(F1_df_Levine_32_ord_size, 
                                            cluster_no = rownames(F1_df_Levine_32_ord_size))
levels(F1_df_Levine_32_ord_size_tidy$cluster_no) <- rownames(F1_df_Levine_32_ord_size)
F1_df_Levine_32_ord_size_tidy <- melt(F1_df_Levine_32_ord_size_tidy, 
                                      id.vars = "cluster_no", 
                                      measure_vars = colnames(F1_df_Levine_32_ord_size), 
                                      variable.name = "method")

F1_df_Levine_13_ord_size <- F1_df_Levine_13[ord_size_Levine_13, ]
F1_df_Levine_13_ord_size_tidy <- data.frame(F1_df_Levine_13_ord_size, 
                                            cluster_no = rownames(F1_df_Levine_13_ord_size))
levels(F1_df_Levine_13_ord_size_tidy$cluster_no) <- rownames(F1_df_Levine_13_ord_size)
F1_df_Levine_13_ord_size_tidy <- melt(F1_df_Levine_13_ord_size_tidy, 
                                      id.vars = "cluster_no", 
                                      measure_vars = colnames(F1_df_Levine_13_ord_size), 
                                      variable.name = "method")

# precision

precision_df_Levine_32_ord_size <- precision_df_Levine_32[ord_size_Levine_32, ]
precision_df_Levine_32_ord_size_tidy <- data.frame(precision_df_Levine_32_ord_size, 
                                                   cluster_no = rownames(precision_df_Levine_32_ord_size))
levels(precision_df_Levine_32_ord_size_tidy$cluster_no) <- rownames(precision_df_Levine_32_ord_size)
precision_df_Levine_32_ord_size_tidy <- melt(precision_df_Levine_32_ord_size_tidy, 
                                             id.vars = "cluster_no", 
                                             measure_vars = colnames(precision_df_Levine_32_ord_size), 
                                             variable.name = "method")

precision_df_Levine_13_ord_size <- precision_df_Levine_13[ord_size_Levine_13, ]
precision_df_Levine_13_ord_size_tidy <- data.frame(precision_df_Levine_13_ord_size, 
                                                   cluster_no = rownames(precision_df_Levine_13_ord_size))
levels(precision_df_Levine_13_ord_size_tidy$cluster_no) <- rownames(precision_df_Levine_13_ord_size)
precision_df_Levine_13_ord_size_tidy <- melt(precision_df_Levine_13_ord_size_tidy, 
                                             id.vars = "cluster_no", 
                                             measure_vars = colnames(precision_df_Levine_13_ord_size), 
                                             variable.name = "method")

# recall

recall_df_Levine_32_ord_size <- recall_df_Levine_32[ord_size_Levine_32, ]
recall_df_Levine_32_ord_size_tidy <- data.frame(recall_df_Levine_32_ord_size, 
                                                cluster_no = rownames(recall_df_Levine_32_ord_size))
levels(recall_df_Levine_32_ord_size_tidy$cluster_no) <- rownames(recall_df_Levine_32_ord_size)
recall_df_Levine_32_ord_size_tidy <- melt(recall_df_Levine_32_ord_size_tidy, 
                                          id.vars = "cluster_no", 
                                          measure_vars = colnames(recall_df_Levine_32_ord_size), 
                                          variable.name = "method")

recall_df_Levine_13_ord_size <- recall_df_Levine_13[ord_size_Levine_13, ]
recall_df_Levine_13_ord_size_tidy <- data.frame(recall_df_Levine_13_ord_size, 
                                                cluster_no = rownames(recall_df_Levine_13_ord_size))
levels(recall_df_Levine_13_ord_size_tidy$cluster_no) <- rownames(recall_df_Levine_13_ord_size)
recall_df_Levine_13_ord_size_tidy <- melt(recall_df_Levine_13_ord_size_tidy, 
                                          id.vars = "cluster_no", 
                                          measure_vars = colnames(recall_df_Levine_13_ord_size), 
                                          variable.name = "method")


# plots of F1 score, precision, and recall for each method, with clusters ordered by size

plot_df_Levine_32_ord_size_tidy <- list(F1_score = F1_df_Levine_32_ord_size_tidy, 
                                        precision = precision_df_Levine_32_ord_size_tidy, 
                                        recall = recall_df_Levine_32_ord_size_tidy)
plot_df_Levine_13_ord_size_tidy <- list(F1_score = F1_df_Levine_13_ord_size_tidy, 
                                        precision = precision_df_Levine_13_ord_size_tidy, 
                                        recall = recall_df_Levine_13_ord_size_tidy)

titles_Levine_32 <- c("F1 score by true cluster size: Levine_2015_marrow_32", 
                      "Precision by true cluster size: Levine_2015_marrow_32", 
                      "Recall by true cluster size: Levine_2015_marrow_32")
titles_Levine_13 <- c("F1 score by true cluster size: Levine_2015_marrow_13", 
                      "Precision by true cluster size: Levine_2015_marrow_13", 
                      "Recall by true cluster size: Levine_2015_marrow_13")

filenames_Levine_32 <- c("../plots/Levine_2015_marrow_32/results_F1_score_by_cluster_size_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_precision_by_cluster_size_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_recall_by_cluster_size_Levine2015marrow32.pdf")
filenames_Levine_13 <- c("../plots/Levine_2015_marrow_13/results_F1_score_by_cluster_size_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_precision_by_cluster_size_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_recall_by_cluster_size_Levine2015marrow13.pdf")


for (i in 1:length(plot_df_Levine_32_ord_size_tidy)) {
  ggplot(plot_df_Levine_32_ord_size_tidy[[i]], 
         aes(x = cluster_no, y = value, group = method, color = method)) + 
    geom_line(stat = "identity") + 
    ylim(0, 1) + 
    xlab("manually gated population number") + 
    ylab(names(plot_df_Levine_32_ord_size_tidy)[i]) + 
    ggtitle(titles_Levine_32[i]) + 
    theme_bw() + 
    theme(legend.key = element_blank())
  
  ggsave(filenames_Levine_32[i], height = 6, width = 8)
}

for (i in 1:length(plot_df_Levine_13_ord_size_tidy)) {
  ggplot(plot_df_Levine_13_ord_size_tidy[[i]], 
         aes(x = cluster_no, y = value, group = method, color = method)) + 
    geom_line(stat = "identity") + 
    ylim(0, 1) + 
    xlab("manually gated population number") + 
    ylab(names(plot_df_Levine_13_ord_size_tidy)[i]) + 
    ggtitle(titles_Levine_13[i]) + 
    theme_bw() + 
    theme(legend.key = element_blank())
  
  ggsave(filenames_Levine_13[i], height = 6, width = 8)
}




#############################################################################
### HEATMAPS: F1 SCORE, PRECISION, AND RECALL SORTED BY TRUE CLUSTER SIZE ###
#############################################################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# rows (true clusters) ranked by size, and columns (methods) arranged in the same order 
# as in main plot of mean F1 scores

df_Levine_32 <- list(F1_df_Levine_32[rev(ord_size_Levine_32), ord_Levine_32], 
                     precision_df_Levine_32[rev(ord_size_Levine_32), ord_Levine_32], 
                     recall_df_Levine_32[rev(ord_size_Levine_32), ord_Levine_32])
df_Levine_13 <- list(F1_df_Levine_13[rev(ord_size_Levine_13), ord_Levine_13], 
                     precision_df_Levine_13[rev(ord_size_Levine_13), ord_Levine_13], 
                     recall_df_Levine_13[rev(ord_size_Levine_13), ord_Levine_13])


titles_Levine_32 <- c("F1 score: Levine_2015_marrow_32", 
                      "Precision: Levine_2015_marrow_32", 
                      "Recall: Levine_2015_marrow_32")
titles_Levine_13 <- c("F1 score: Levine_2015_marrow_13", 
                      "Precision: Levine_2015_marrow_13", 
                      "Recall: Levine_2015_marrow_13")

filenames_Levine_32 <- c("../plots/Levine_2015_marrow_32/results_heatmap_F1_score_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_heatmap_precision_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_heatmap_recall_Levine2015marrow32.pdf")
filenames_Levine_13 <- c("../plots/Levine_2015_marrow_13/results_heatmap_F1_score_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_heatmap_precision_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_heatmap_recall_Levine2015marrow13.pdf")


# heatmaps of F1 score, precision, and recall for each true cluster

for (i in 1:length(df_Levine_32)) {
  pheatmap(df_Levine_32[[i]], color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
           breaks = seq(0, 1, length.out = 100), display_numbers = TRUE, 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           main = titles_Levine_32[i], filename = filenames_Levine_32[i], 
           width = 6, height = 6)
}

for (i in 1:length(df_Levine_13)) {
  pheatmap(df_Levine_13[[i]], color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
           breaks = seq(0, 1, length.out = 100), display_numbers = TRUE, 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           main = titles_Levine_13[i], filename = filenames_Levine_13[i], 
           width = 6, height = 9)
}




##########################################################################
### DETECTION OF RARE CELL POPULATION: F1 SCORE, PRECISION, AND RECALL ###
##########################################################################

### for data sets with a single rare cell population of interest (Mosmann_2014_rare)


# arrange by F1 score

ord_Mosmann <- rev(order(F1_df_Mosmann))

precision_df_Mosmann_ord <- unlist(precision_df_Mosmann[ord_Mosmann])
recall_df_Mosmann_ord <- unlist(recall_df_Mosmann[ord_Mosmann])
F1_df_Mosmann_ord <- unlist(F1_df_Mosmann[ord_Mosmann])


# tidy data format (for ggplot)

plot_data_Mosmann <- data.frame(precision = precision_df_Mosmann_ord, 
                                recall = recall_df_Mosmann_ord, 
                                F1_score = F1_df_Mosmann_ord)
plot_data_Mosmann["method"] <- factor(rownames(plot_data_Mosmann), 
                                      levels = rownames(plot_data_Mosmann))
plot_data_Mosmann <- melt(plot_data_Mosmann, 
                          id.vars = "method", 
                          measure.vars = c("F1_score", "precision", "recall"))

# bar plots

ggplot(plot_data_Mosmann, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Rare cell population: Mosmann_2014_rare") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Mosmann_2014_rare/results_F1_precision_recall_Mosmann2014rare.pdf", height = 6, width = 7)


