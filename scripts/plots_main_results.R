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

res_Nilsson <- list(ACCENSE = res_ACCENSE_Nilsson, 
                    DensVM = res_DensVM_Nilsson, 
                    FLOCK = res_FLOCK_Nilsson, 
                    flowMeans = res_flowMeans_Nilsson, 
                    FlowSOM = res_FlowSOM_Nilsson, 
                    FlowSOM_meta = res_FlowSOM_meta_Nilsson, 
                    immunoClust = res_immunoClust_Nilsson, 
                    immunoClust_all = res_immunoClust_all_Nilsson, 
                    kmeans = res_kmeans_Nilsson, 
                    PhenoGraph = res_PhenoGraph_Nilsson, 
                    Rclusterpp = res_Rclusterpp_Nilsson, 
                    SamSPECTRAL = res_SamSPECTRAL_Nilsson, 
                    SWIFT = res_SWIFT_Nilsson)

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
precision_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$pr))
precision_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$pr))

recall_df_Levine_32 <- sapply(res_Levine_32, function(res) res$re)
recall_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$re))
recall_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$re))
recall_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$re))

F1_df_Levine_32 <- sapply(res_Levine_32, function(res) res$F1)
F1_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$F1))
F1_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$F1))
F1_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$F1))

labels_matched_df_Levine_32 <- sapply(res_Levine_32, function(res) res$labels_matched)
labels_matched_df_Levine_13 <- collapse_data_frame_NAs(sapply(res_Levine_13, function(res) res$labels_matched))  # NA = no match
labels_matched_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$labels_matched))
labels_matched_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$labels_matched))

n_cells_df_Levine_32 <- sapply(res_Levine_32, function(res) res$n_cells)
n_cells_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$n_cells))
n_cells_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$n_cells))
n_cells_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$n_cells))




#################################
### MEAN F1 SCORE: UNWEIGHTED ###
#################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# number of methods

n_methods_Levine_32 <- ncol(n_cells_df_Levine_32)
n_methods_Levine_13 <- ncol(n_cells_df_Levine_13)

# mean F1 score across all true clusters (unweighted)

mean_F1_Levine_32 <- colMeans(F1_df_Levine_32)
mean_F1_Levine_13 <- colMeans(F1_df_Levine_13)

# arrange in descending order

ord_Levine_32 <- rev(order(mean_F1_Levine_32))
ord_Levine_13 <- rev(order(mean_F1_Levine_13))

mean_F1_Levine_32_ord <- mean_F1_Levine_32[ord_Levine_32]
mean_F1_Levine_13_ord <- mean_F1_Levine_13[ord_Levine_13]

mean_F1_Levine_32_ord
mean_F1_Levine_13_ord


# tidy data format (for ggplot)

mean_F1_Levine_32_tidy <- data.frame(value = mean_F1_Levine_32_ord)
mean_F1_Levine_32_tidy["method"] <- factor(rownames(mean_F1_Levine_32_tidy), 
                                           levels = rownames(mean_F1_Levine_32_tidy))
mean_F1_Levine_32_tidy

mean_F1_Levine_13_tidy <- data.frame(value = mean_F1_Levine_13_ord)
mean_F1_Levine_13_tidy["method"] <- factor(rownames(mean_F1_Levine_13_tidy), 
                                           levels = rownames(mean_F1_Levine_13_tidy))
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




###################################################################
### MEAN F1 SCORE: WEIGHTED BY NUMBER OF CELLS PER TRUE CLUSTER ###
###################################################################

### continued from above


# number of cells per true clusters

n_cells_truth_Levine_32 <- matrix(rep(tbl_truth_Levine_32, n_methods_Levine_32), ncol = n_methods_Levine_32)
n_cells_truth_Levine_13 <- matrix(rep(tbl_truth_Levine_13, n_methods_Levine_13), ncol = n_methods_Levine_13)

n_cells_truth_Levine_32
n_cells_truth_Levine_13

# mean F1 score across all true clusters (weighted by number of cells in true cluster)

mean_F1_Levine_32_weighted <- colSums(F1_df_Levine_32 * n_cells_truth_Levine_32) / colSums(n_cells_truth_Levine_32)
mean_F1_Levine_13_weighted <- colSums(F1_df_Levine_13 * n_cells_truth_Levine_13) / colSums(n_cells_truth_Levine_13)

# arrange in descending order

ord_Levine_32_weighted <- rev(order(mean_F1_Levine_32_weighted))
ord_Levine_13_weighted <- rev(order(mean_F1_Levine_13_weighted))

mean_F1_Levine_32_ord_weighted <- mean_F1_Levine_32_weighted[ord_Levine_32_weighted]
mean_F1_Levine_13_ord_weighted <- mean_F1_Levine_13_weighted[ord_Levine_13_weighted]

mean_F1_Levine_32_ord_weighted
mean_F1_Levine_13_ord_weighted


# tidy data format (for ggplot)

mean_F1_Levine_32_tidy_weighted <- data.frame(value = mean_F1_Levine_32_ord_weighted)
mean_F1_Levine_32_tidy_weighted["method"] <- factor(rownames(mean_F1_Levine_32_tidy_weighted), 
                                                    levels = rownames(mean_F1_Levine_32_tidy_weighted))
mean_F1_Levine_32_tidy_weighted

mean_F1_Levine_13_tidy_weighted <- data.frame(value = mean_F1_Levine_13_ord_weighted)
mean_F1_Levine_13_tidy_weighted["method"] <- factor(rownames(mean_F1_Levine_13_tidy_weighted), 
                                                    levels = rownames(mean_F1_Levine_13_tidy_weighted))
mean_F1_Levine_13_tidy_weighted


# bar plots

ggplot(mean_F1_Levine_32_tidy_weighted, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.2f", round(value, 2)), y = value + 0.03), size = 3.5) + 
  ylim(0, 1) + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score (weighted): Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/Levine_2015_marrow_32/results_mean_F1_score_weighted_Levine2015marrow32.pdf", height = 6, width = 5.5)


ggplot(mean_F1_Levine_13_tidy_weighted, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.2f", round(value, 2)), y = value + 0.03), size = 3.5) + 
  ylim(0, 1) + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score (weighted): Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/Levine_2015_marrow_13/results_mean_F1_score_weighted_Levine2015marrow13.pdf", height = 6, width = 5.5)




############################################
### MEAN F1 SCORE, PRECISION, AND RECALL ###
############################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# mean precision and mean recall across all true clusters (unweighted)

mean_precision_Levine_32 <- colMeans(precision_df_Levine_32)
mean_precision_Levine_13 <- colMeans(precision_df_Levine_13)

mean_recall_Levine_32 <- colMeans(recall_df_Levine_32)
mean_recall_Levine_13 <- colMeans(recall_df_Levine_13)

# arrange in descending order of mean F1 score (unweighted)

mean_precision_Levine_32_ord <- mean_precision_Levine_32[ord_Levine_32]
mean_precision_Levine_13_ord <- mean_precision_Levine_13[ord_Levine_13]

mean_recall_Levine_32_ord <- mean_recall_Levine_32[ord_Levine_32]
mean_recall_Levine_13_ord <- mean_recall_Levine_13[ord_Levine_13]

mean_precision_Levine_32_ord
mean_precision_Levine_13_ord

mean_recall_Levine_32_ord
mean_recall_Levine_13_ord


# tidy data format (for ggplot)

plot_data_Levine_32 <- data.frame(F1_score = mean_F1_Levine_32_ord, 
                                  precision = mean_precision_Levine_32_ord, 
                                  recall = mean_recall_Levine_32_ord)
plot_data_Levine_32["method"] <- factor(rownames(plot_data_Levine_32), 
                                        levels = rownames(plot_data_Levine_32))
plot_data_Levine_32 <- melt(plot_data_Levine_32, 
                            id.vars = "method", 
                            measure.vars = c("F1_score", "precision", "recall"))

plot_data_Levine_13 <- data.frame(F1_score = mean_F1_Levine_13_ord, 
                                  precision = mean_precision_Levine_13_ord, 
                                  recall = mean_recall_Levine_13_ord)
plot_data_Levine_13["method"] <- factor(rownames(plot_data_Levine_13), 
                                        levels = rownames(plot_data_Levine_13))
plot_data_Levine_13 <- melt(plot_data_Levine_13, 
                            id.vars = "method", 
                            measure.vars = c("F1_score", "precision", "recall"))


# bar plots of mean F1 score, mean precision, and mean recall (in same order as previously)

ggplot(plot_data_Levine_32, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Mean F1 score, precision, and recall: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Levine_2015_marrow_32/results_mean_F1_precision_recall_Levine2015marrow32.pdf", height = 6, width = 7)


ggplot(plot_data_Levine_13, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Mean F1 score, precision, and recall: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Levine_2015_marrow_13/results_mean_F1_precision_recall_Levine2015marrow13.pdf", height = 6, width = 7)



##########################################################################################
### MEAN F1 SCORE, PRECISION, AND RECALL: WEIGHTED BY NUMBER OF CELLS PER TRUE CLUSTER ###
##########################################################################################

### continued from above


# mean precision and mean recall across all true clusters (weighted by number of cells in true cluster)

mean_precision_Levine_32_weighted <- colSums(precision_df_Levine_32 * n_cells_truth_Levine_32) / colSums(n_cells_truth_Levine_32)
mean_precision_Levine_13_weighted <- colSums(precision_df_Levine_13 * n_cells_truth_Levine_13) / colSums(n_cells_truth_Levine_13)

mean_recall_Levine_32_weighted <- colSums(recall_df_Levine_32 * n_cells_truth_Levine_32) / colSums(n_cells_truth_Levine_32)
mean_recall_Levine_13_weighted <- colSums(recall_df_Levine_13 * n_cells_truth_Levine_13) / colSums(n_cells_truth_Levine_13)

# arrange in descending order of mean F1 score (weighted by number of cells in true cluster)

mean_precision_Levine_32_ord_weighted <- mean_precision_Levine_32_weighted[ord_Levine_32_weighted]
mean_precision_Levine_13_ord_weighted <- mean_precision_Levine_13_weighted[ord_Levine_13_weighted]

mean_recall_Levine_32_ord_weighted <- mean_recall_Levine_32_weighted[ord_Levine_32_weighted]
mean_recall_Levine_13_ord_weighted <- mean_recall_Levine_13_weighted[ord_Levine_13_weighted]

mean_precision_Levine_32_ord_weighted
mean_precision_Levine_13_ord_weighted

mean_recall_Levine_32_ord_weighted
mean_recall_Levine_13_ord_weighted


# tidy data format (for ggplot)

plot_data_Levine_32_weighted <- data.frame(F1_score = mean_F1_Levine_32_ord_weighted, 
                                           precision = mean_precision_Levine_32_ord_weighted, 
                                           recall = mean_recall_Levine_32_ord_weighted)
plot_data_Levine_32_weighted["method"] <- factor(rownames(plot_data_Levine_32_weighted), 
                                                 levels = rownames(plot_data_Levine_32_weighted))
plot_data_Levine_32_weighted <- melt(plot_data_Levine_32_weighted, 
                                     id.vars = "method", 
                                     measure.vars = c("F1_score", "precision", "recall"))

plot_data_Levine_13_weighted <- data.frame(F1_score = mean_F1_Levine_13_ord_weighted, 
                                           precision = mean_precision_Levine_13_ord_weighted, 
                                           recall = mean_recall_Levine_13_ord_weighted)
plot_data_Levine_13_weighted["method"] <- factor(rownames(plot_data_Levine_13_weighted), 
                                                 levels = rownames(plot_data_Levine_13_weighted))
plot_data_Levine_13_weighted <- melt(plot_data_Levine_13_weighted, 
                                     id.vars = "method", 
                                     measure.vars = c("F1_score", "precision", "recall"))


# bar plots of mean F1 score, mean precision, and mean recall (in same order as previously)

ggplot(plot_data_Levine_32_weighted, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Mean F1 score, precision, and recall (weighted): Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Levine_2015_marrow_32/results_mean_F1_precision_recall_weighted_Levine2015marrow32.pdf", height = 6, width = 7)


ggplot(plot_data_Levine_13_weighted, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Mean F1 score, precision, and recall (weighted): Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Levine_2015_marrow_13/results_mean_F1_precision_recall_weighted_Levine2015marrow13.pdf", height = 6, width = 7)




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




#############################################################################
### HEATMAPS: F1 SCORE, PRECISION, AND RECALL SORTED BY TRUE CLUSTER SIZE ###
#############################################################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# rows (true clusters) ranked by size, and columns (methods) arranged in the same order 
# as in main plot of mean F1 scores

ord_size_rev_Levine_32 <- rev(order(tbl_truth_Levine_32))
ord_size_rev_Levine_13 <- rev(order(tbl_truth_Levine_13))


df_Levine_32 <- list(F1_df_Levine_32[ord_size_rev_Levine_32, ord_Levine_32], 
                     precision_df_Levine_32[ord_size_rev_Levine_32, ord_Levine_32], 
                     recall_df_Levine_32[ord_size_rev_Levine_32, ord_Levine_32])

df_Levine_13 <- list(F1_df_Levine_13[ord_size_rev_Levine_13, ord_Levine_13], 
                     precision_df_Levine_13[ord_size_rev_Levine_13, ord_Levine_13], 
                     recall_df_Levine_13[ord_size_rev_Levine_13, ord_Levine_13])


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

### for data sets with a single rare cell population of interest (Nilsson_2013_HSC, Mosmann_2014_activ)


# arrange by F1 score

ord_Nilsson <- rev(order(F1_df_Nilsson))
ord_Mosmann <- rev(order(F1_df_Mosmann))


precision_df_Nilsson_ord <- unlist(precision_df_Nilsson[ord_Nilsson])
precision_df_Mosmann_ord <- unlist(precision_df_Mosmann[ord_Mosmann])

recall_df_Nilsson_ord <- unlist(recall_df_Nilsson[ord_Nilsson])
recall_df_Mosmann_ord <- unlist(recall_df_Mosmann[ord_Mosmann])

F1_df_Nilsson_ord <- unlist(F1_df_Nilsson[ord_Nilsson])
F1_df_Mosmann_ord <- unlist(F1_df_Mosmann[ord_Mosmann])


# tidy data format (for ggplot)

plot_data_Nilsson <- data.frame(precision = precision_df_Nilsson_ord, 
                                recall = recall_df_Nilsson_ord, 
                                F1_score = F1_df_Nilsson_ord)
plot_data_Nilsson["method"] <- factor(rownames(plot_data_Nilsson), 
                                      levels = rownames(plot_data_Nilsson))
plot_data_Nilsson <- melt(plot_data_Nilsson, 
                          id.vars = "method", 
                          measure.vars = c("F1_score", "precision", "recall"))

plot_data_Mosmann <- data.frame(precision = precision_df_Mosmann_ord, 
                                recall = recall_df_Mosmann_ord, 
                                F1_score = F1_df_Mosmann_ord)
plot_data_Mosmann["method"] <- factor(rownames(plot_data_Mosmann), 
                                      levels = rownames(plot_data_Mosmann))
plot_data_Mosmann <- melt(plot_data_Mosmann, 
                          id.vars = "method", 
                          measure.vars = c("F1_score", "precision", "recall"))

# bar plots

ggplot(plot_data_Nilsson, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Rare cell population: Nilsson_2013_HSC") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Nilsson_2013_HSC/results_F1_precision_recall_Nilsson2013HSC.pdf", height = 6, width = 7)


ggplot(plot_data_Mosmann, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  ggtitle("Rare cell population: Mosmann_2014_activ") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/Mosmann_2014_activ/results_F1_precision_recall_Mosmann2014activ.pdf", height = 6, width = 7)


