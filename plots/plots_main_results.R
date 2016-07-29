#########################################################################################
# R script to generate main plots of results
#
# Lukas M. Weber, March 2016
#########################################################################################


library(ggplot2)
library(reshape2)
library(cowplot)  # note masks ggplot2::ggsave
library(ggrepel)

# helper function
source("helper_collapse_data_frame.R")

# load results directories (depending on whether automatic or manually selected number of clusters)
source("load_results_directories.R")

# load results from previous steps
source("load_results_ACCENSE.R")
source("load_results_DensVM.R")
source("load_results_FLOCK.R")
source("load_results_PhenoGraph.R")
source("load_results_Rclusterpp.R")
source("load_results_SWIFT.R")
source("load_results_truth.R")
source("load_results_all_other_methods.R")

# load runtime results
source("load_results_runtime.R")

# color-blind friendly palettes (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cb_pal_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_pal_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ggplot2 default palette (access with ggplot_build(p)$data)
gg_pal <- c("#F8766D", "#00BA38", "#619CFF")




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
                    Rclusterpp = res_Rclusterpp_Mosmann, 
                    SamSPECTRAL = res_SamSPECTRAL_Mosmann, 
                    SWIFT = res_SWIFT_Mosmann)


# collapse into data frames (using helper functions to pad missing values with zeros or NAs)

precision_df_Levine_32 <- as.data.frame(sapply(res_Levine_32, function(res) res$pr))
precision_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$pr))
precision_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$pr))
precision_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$pr))

recall_df_Levine_32 <- as.data.frame(sapply(res_Levine_32, function(res) res$re))
recall_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$re))
recall_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$re))
recall_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$re))

F1_df_Levine_32 <- as.data.frame(sapply(res_Levine_32, function(res) res$F1))
F1_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$F1))
F1_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$F1))
F1_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$F1))

labels_matched_df_Levine_32 <- as.data.frame(sapply(res_Levine_32, function(res) res$labels_matched))
labels_matched_df_Levine_13 <- collapse_data_frame_NAs(sapply(res_Levine_13, function(res) res$labels_matched))  # NA = no match
labels_matched_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$labels_matched))
labels_matched_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$labels_matched))

n_cells_df_Levine_32 <- as.data.frame(sapply(res_Levine_32, function(res) res$n_cells))
n_cells_df_Levine_13 <- collapse_data_frame_zeros(sapply(res_Levine_13, function(res) res$n_cells))
n_cells_df_Nilsson <- as.data.frame(lapply(res_Nilsson, function(res) res$n_cells))
n_cells_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$n_cells))




#####################
### MEAN F1 SCORE ###
#####################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


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

barplot_mean_F1_Levine_32 <- 
  ggplot(mean_F1_Levine_32_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "royalblue3") + 
  geom_text(aes(label = sprintf("%.3f", round(value, 3)), y = value + 0.08, angle = 90), size = 3.5) + 
  ylim(0, 1) + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

barplot_mean_F1_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/results_barplot_mean_F1_Levine2015marrow32.pdf", 
                width = 5, height = 5)


barplot_mean_F1_Levine_13 <- 
  ggplot(mean_F1_Levine_13_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "royalblue3") + 
  geom_text(aes(label = sprintf("%.3f", round(value, 3)), y = value + 0.08, angle = 90), size = 3.5) + 
  ylim(0, 1) + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

barplot_mean_F1_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/results_barplot_mean_F1_Levine2015marrow13.pdf", 
                width = 5, height = 5)




###########################
### F1 SCORE: BOX PLOTS ###
###########################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# arrange in same order as previously

F1_df_Levine_32_ord <- F1_df_Levine_32[, ord_Levine_32]
F1_df_Levine_13_ord <- F1_df_Levine_13[, ord_Levine_13]

F1_df_Levine_32_ord
F1_df_Levine_13_ord


# tidy data format (for ggplot)

F1_df_Levine_32_tidy <- data.frame(value = as.vector(as.matrix(F1_df_Levine_32_ord)))
F1_df_Levine_32_tidy["method"] <- rep(factor(colnames(F1_df_Levine_32_ord), 
                                             levels = colnames(F1_df_Levine_32_ord)), 
                                      each = nrow(F1_df_Levine_32_ord))
F1_df_Levine_32_tidy

F1_df_Levine_13_tidy <- data.frame(value = as.vector(as.matrix(F1_df_Levine_13_ord)))
F1_df_Levine_13_tidy["method"] <- rep(factor(colnames(F1_df_Levine_13_ord), 
                                             levels = colnames(F1_df_Levine_13_ord)), 
                                      each = nrow(F1_df_Levine_13_ord))
F1_df_Levine_13_tidy


# box plots

boxplots_F1_Levine_32 <- 
  ggplot(F1_df_Levine_32_tidy, aes(x = method, y = value)) + 
  geom_boxplot(col = "gray50", fill = "aliceblue") + 
  geom_point(shape = 1, col = "darkblue") + 
  stat_summary(fun.y = mean, color = "red", geom = "point", shape = 1) + 
  ylim(0, 1) + 
  ylab("F1 score") + 
  ggtitle("F1 score: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

boxplots_F1_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/results_boxplots_F1_Levine2015marrow32.pdf", 
                width = 5, height = 5)


boxplots_F1_Levine_13 <- 
  ggplot(F1_df_Levine_13_tidy, aes(x = method, y = value)) + 
  geom_boxplot(col = "gray50", fill = "aliceblue") + 
  geom_point(shape = 1, col = "darkblue") + 
  stat_summary(fun.y = mean, color = "red", geom = "point", shape = 1) + 
  ylim(0, 1) + 
  ylab("F1 score") + 
  ggtitle("F1 score: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

boxplots_F1_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/results_boxplots_F1_Levine2015marrow13.pdf", 
                width = 5, height = 5)




############################################
### MEAN F1 SCORE, PRECISION, AND RECALL ###
############################################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


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

barplot_mean_F1_pr_re_Levine_32 <- 
  ggplot(plot_data_Levine_32, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
  ylim(0, 1) + 
  ylab("") + 
  ggtitle("Mean F1 score, precision, and recall: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.73, 0.95), 
        legend.direction = "horizontal", 
        legend.key.size = unit(4, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

barplot_mean_F1_pr_re_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/results_barplot_mean_F1_pr_re_Levine2015marrow32.pdf", 
                width = 5, height = 5)


barplot_mean_F1_pr_re_Levine_13 <- 
  ggplot(plot_data_Levine_13, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
  ylim(0, 1) + 
  ylab("") + 
  ggtitle("Mean F1 score, precision, and recall: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.73, 0.95), 
        legend.direction = "horizontal", 
        legend.key.size = unit(4, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

barplot_mean_F1_pr_re_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/results_barplot_mean_F1_pr_re_Levine2015marrow13.pdf", 
                width = 5, height = 5)




########################
### POPULATION SIZES ###
########################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# tidy data format (for ggplot)

n_cells_truth_Levine_32_tidy <- as.data.frame(tbl_truth_Levine_32)
n_cells_truth_Levine_13_tidy <- as.data.frame(tbl_truth_Levine_13)

colnames(n_cells_truth_Levine_32_tidy) <- c("population", "value")
colnames(n_cells_truth_Levine_13_tidy) <- c("population", "value")


# plot number of cells in each true (manually gated) population

plot_n_cells_Levine_32 <- 
  ggplot(n_cells_truth_Levine_32_tidy, aes(x = population, y = value)) + 
  geom_bar(stat = "identity", fill = "darkgray") + 
  geom_text(aes(label = value, y = value + 600, angle = 90), hjust = "left", size = 3.5) + 
  ylim(0, 30000) + 
  xlab("manually gated population") + 
  ylab("number of cells") + 
  ggtitle("Manually gated populations: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12))

plot_n_cells_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/results_no_of_cells_Levine2015marrow32.pdf", 
                width = 5, height = 5)


plot_n_cells_Levine_13 <- 
  ggplot(n_cells_truth_Levine_13_tidy, aes(x = population, y = value)) + 
  geom_bar(stat = "identity", fill = "darkgray") + 
  geom_text(aes(label = value, y = value + 300, angle = 90), hjust = "left", size = 3.5) + 
  ylim(0, 15750) + 
  xlab("manually gated population") + 
  ylab("number of cells") + 
  ggtitle("Manually gated populations: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.text.x = element_text(size = 9))

plot_n_cells_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/results_no_of_cells_Levine2015marrow13.pdf", 
                width = 5, height = 5)




##########################################################################
### DETECTION OF RARE CELL POPULATION: F1 SCORE, PRECISION, AND RECALL ###
##########################################################################

# for data sets with a single rare cell population of interest (Nilsson_2013_HSC, Mosmann_2014_activ)


# arrange by F1 score

ord_Nilsson <- rev(order(F1_df_Nilsson))
ord_Mosmann <- rev(order(F1_df_Mosmann))


precision_df_Nilsson_ord <- unlist(precision_df_Nilsson[ord_Nilsson])
precision_df_Mosmann_ord <- unlist(precision_df_Mosmann[ord_Mosmann])

recall_df_Nilsson_ord <- unlist(recall_df_Nilsson[ord_Nilsson])
recall_df_Mosmann_ord <- unlist(recall_df_Mosmann[ord_Mosmann])

F1_df_Nilsson_ord <- unlist(F1_df_Nilsson[ord_Nilsson])
F1_df_Mosmann_ord <- unlist(F1_df_Mosmann[ord_Mosmann])

F1_df_Nilsson_ord
F1_df_Mosmann_ord


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

barplot_F1_pr_re_Nilsson <- 
  ggplot(plot_data_Nilsson, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
  ylim(0, 1.05) + 
  ylab("") + 
  ggtitle("Rare cell population: Nilsson_2013_HSC") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.73, 0.955), 
        legend.direction = "horizontal", 
        legend.key.size = unit(4, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

barplot_F1_pr_re_Nilsson
ggplot2::ggsave("../plots/Nilsson_2013_HSC/results_barplot_F1_pr_re_Nilsson2013HSC.pdf", 
                width = 5, height = 5)


barplot_F1_pr_re_Mosmann <- 
  ggplot(plot_data_Mosmann, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
  ylim(0, 1.05) + 
  ylab("") + 
  ggtitle("Rare cell population: Mosmann_2014_activ") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.73, 0.955), 
        legend.direction = "horizontal", 
        legend.key.size = unit(4, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

barplot_F1_pr_re_Mosmann
ggplot2::ggsave("../plots/Mosmann_2014_activ/results_barplot_F1_pr_re_Mosmann2014activ.pdf", 
                width = 5, height = 5)




##########################
### RUNTIME: BAR PLOTS ###
##########################

# runtime results are loaded with a separate script (sourced at the beginning of this script)


# convert data frames to tidy data format (for ggplot)

runtime_Levine_32_tidy <- data.frame(value = runtime_Levine_32_ord)
runtime_Levine_32_tidy["method"] <- factor(rownames(runtime_Levine_32_tidy), 
                                           levels = rownames(runtime_Levine_32_tidy))

runtime_Levine_13_tidy <- data.frame(value = runtime_Levine_13_ord)
runtime_Levine_13_tidy["method"] <- factor(rownames(runtime_Levine_13_tidy), 
                                           levels = rownames(runtime_Levine_13_tidy))

runtime_Nilsson_tidy <- data.frame(value = runtime_Nilsson_ord)
runtime_Nilsson_tidy["method"] <- factor(rownames(runtime_Nilsson_tidy), 
                                         levels = rownames(runtime_Nilsson_tidy))

runtime_Mosmann_tidy <- data.frame(value = runtime_Mosmann_ord)
runtime_Mosmann_tidy["method"] <- factor(rownames(runtime_Mosmann_tidy), 
                                         levels = rownames(runtime_Mosmann_tidy))


# single or multiple cores used

runtime_Levine_32_tidy$cores <- "single core"
runtime_Levine_32_tidy[c("SWIFT", "Rclusterpp"), "cores"] <- "multiple cores"
runtime_Levine_32_tidy$cores <- factor(runtime_Levine_32_tidy$cores, 
                                       levels = c("single core", "multiple cores"))

runtime_Levine_13_tidy$cores <- "single core"
runtime_Levine_13_tidy[c("SWIFT", "Rclusterpp"), "cores"] <- "multiple cores"
runtime_Levine_13_tidy$cores <- factor(runtime_Levine_13_tidy$cores, 
                                       levels = c("single core", "multiple cores"))

runtime_Nilsson_tidy$cores <- "single core"
runtime_Nilsson_tidy[c("SWIFT", "Rclusterpp"), "cores"] <- "multiple cores"
runtime_Nilsson_tidy$cores <- factor(runtime_Nilsson_tidy$cores, 
                                     levels = c("single core", "multiple cores"))

runtime_Mosmann_tidy$cores <- "single core"
runtime_Mosmann_tidy[c("SWIFT", "Rclusterpp"), "cores"] <- "multiple cores"
runtime_Mosmann_tidy$cores <- factor(runtime_Mosmann_tidy$cores, 
                                     levels = c("single core", "multiple cores"))


# check

runtime_Levine_32_tidy
runtime_Levine_13_tidy
runtime_Nilsson_tidy
runtime_Mosmann_tidy


# bar plots

runtime_barplot_Levine_32 <- 
  ggplot(runtime_Levine_32_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", aes(fill = cores)) + 
  scale_fill_manual(values = c("mediumpurple", "darkgray")) + 
  geom_text(aes(label = round(value, 0), y = value + 800, angle = 90), hjust = "left", size = 3.5) + 
  ggtitle("Runtime: Levine_2015_marrow_32") + 
  ylim(0, 38000) + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.16, 0.91), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

runtime_barplot_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/runtime_barplot_Levine2015marrow32.pdf", 
                width = 5, height = 5)


runtime_barplot_Levine_13 <- 
  ggplot(runtime_Levine_13_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", aes(fill = cores)) + 
  scale_fill_manual(values = c("mediumpurple", "darkgray")) + 
  geom_text(aes(label = round(value, 0), y = value + 200, angle = 90), hjust = "left", size = 3.5) + 
  ggtitle("Runtime: Levine_2015_marrow_13") + 
  ylim(0, 10500) + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.16, 0.91), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

runtime_barplot_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/runtime_barplot_Levine2015marrow13.pdf", 
                width = 5, height = 5)


runtime_barplot_Nilsson <- 
  ggplot(runtime_Nilsson_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", aes(fill = cores)) + 
  scale_fill_manual(values = c("mediumpurple", "darkgray")) + 
  geom_text(aes(label = round(value, 0), y = value + 200, angle = 90), hjust = "left", size = 3.5) + 
  ggtitle("Runtime: Nilsson_2013_HSC") + 
  ylim(0, 10000) + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.16, 0.91), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

runtime_barplot_Nilsson
ggplot2::ggsave("../plots/Nilsson_2013_HSC/runtime_barplot_Nilsson2013HSC.pdf", 
                width = 5, height = 5)


runtime_barplot_Mosmann <- 
  ggplot(runtime_Mosmann_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", aes(fill = cores)) + 
  scale_fill_manual(values = c("mediumpurple", "darkgray")) + 
  geom_text(aes(label = round(value, 0), y = value + 350, angle = 90), hjust = "left", size = 3.5) + 
  ggtitle("Runtime: Mosmann_2014_activ") + 
  ylim(0, 16750) + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.16, 0.91), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

runtime_barplot_Mosmann
ggplot2::ggsave("../plots/Mosmann_2014_activ/runtime_barplot_Mosmann2014activ.pdf", 
                width = 5, height = 5)




##########################################################
### RUNTIME: SCATTER PLOTS (RUNTIME VS. MEAN F1 SCORE) ###
##########################################################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# create tidy data frames

runtime_vs_F1_Levine_32_tidy <- runtime_Levine_32_tidy
colnames(runtime_vs_F1_Levine_32_tidy)[1] <- "runtime"
runtime_vs_F1_Levine_32_tidy$mean_F1 <- mean_F1_Levine_32[rownames(runtime_vs_F1_Levine_32_tidy)]

runtime_vs_F1_Levine_13_tidy <- runtime_Levine_13_tidy
colnames(runtime_vs_F1_Levine_13_tidy)[1] <- "runtime"
runtime_vs_F1_Levine_13_tidy$mean_F1 <- mean_F1_Levine_13[rownames(runtime_vs_F1_Levine_13_tidy)]

# check

runtime_vs_F1_Levine_32_tidy
runtime_vs_F1_Levine_13_tidy


# scatter plots

runtime_scatterplot_Levine_32 <- 
  ggplot(runtime_vs_F1_Levine_32_tidy, aes(x = mean_F1, y = runtime)) + 
  geom_point(shape = 4, size = 2, stroke = 1, color = "darkorchid4") + 
  geom_text_repel(aes(label = method), size = 2.5, box.padding = unit(0.3, "lines")) + 
  xlim(0.15, 0.8) + 
  ylim(-1000, 38000) + 
  ggtitle("Runtime vs. mean F1: Levine_2015_marrow_32") + 
  xlab("mean F1 score") + 
  ylab("runtime (seconds)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12))

runtime_scatterplot_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/runtime_scatterplot_Levine2015marrow32.pdf", 
                width = 5, height = 5)


runtime_scatterplot_Levine_13 <- 
  ggplot(runtime_vs_F1_Levine_13_tidy, aes(x = mean_F1, y = runtime)) + 
  geom_point(shape = 4, size = 2, stroke = 1, color = "darkorchid4") + 
  geom_text_repel(aes(label = method), size = 2.5, box.padding = unit(0.3, "lines")) + 
  xlim(0.15, 0.8) + 
  ylim(-1000, 10500) + 
  ggtitle("Runtime vs. mean F1: Levine_2015_marrow_13") + 
  xlab("mean F1 score") + 
  ylab("runtime (seconds)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12))

runtime_scatterplot_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/runtime_scatterplot_Levine2015marrow13.pdf", 
                width = 5, height = 5)


# for data sets with a single rare cell population of interest (Nilsson_2013_HSC, Mosmann_2014_activ)

runtime_vs_F1_Nilsson_tidy <- runtime_Nilsson_tidy
colnames(runtime_vs_F1_Nilsson_tidy)[1] <- "runtime"
runtime_vs_F1_Nilsson_tidy$F1 <- F1_df_Nilsson_ord[rownames(runtime_vs_F1_Nilsson_tidy)]

runtime_vs_F1_Mosmann_tidy <- runtime_Mosmann_tidy
colnames(runtime_vs_F1_Mosmann_tidy)[1] <- "runtime"
runtime_vs_F1_Mosmann_tidy$F1 <- F1_df_Mosmann_ord[rownames(runtime_vs_F1_Mosmann_tidy)]

# check

runtime_vs_F1_Nilsson_tidy
runtime_vs_F1_Mosmann_tidy


# scatter plots

runtime_scatterplot_Nilsson <- 
  ggplot(runtime_vs_F1_Nilsson_tidy, aes(x = F1, y = runtime)) + 
  geom_point(shape = 4, size = 2, stroke = 1, color = "darkorchid4") + 
  geom_text_repel(aes(label = method), size = 2.5, box.padding = unit(0.3, "lines"), force = 3) + 
  xlim(-0.1, 0.75) + 
  ylim(-1000, 10000) + 
  ggtitle("Runtime vs. F1: Nilsson_2013_HSC") + 
  xlab("F1 score") + 
  ylab("runtime (seconds)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12))

runtime_scatterplot_Nilsson
ggplot2::ggsave("../plots/Nilsson_2013_HSC/runtime_scatterplot_Nilsson2013HSC.pdf", 
                width = 5, height = 5)


runtime_scatterplot_Mosmann <- 
  ggplot(runtime_vs_F1_Mosmann_tidy, aes(x = F1, y = runtime)) + 
  geom_point(shape = 4, size = 2, stroke = 1, color = "darkorchid4") + 
  geom_text_repel(aes(label = method), size = 2.5, box.padding = unit(0.45, "lines")) + 
  xlim(-0.1, 0.75) + 
  ylim(-1000, 16750) + 
  ggtitle("Runtime vs. F1: Mosmann_2014_activ") + 
  xlab("F1 score") + 
  ylab("runtime (seconds)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12))

runtime_scatterplot_Mosmann
ggplot2::ggsave("../plots/Mosmann_2014_activ/runtime_scatterplot_Mosmann2014activ.pdf", 
                width = 5, height = 5)




#########################
### MULTI-PANEL PLOTS ###
#########################

# combine into one multi-panel plot for each data set


# Levine_2015_marrow_32
ggdraw() + 
  draw_plot(barplot_mean_F1_Levine_32, 0.05, 0.66, 0.4, 0.33) + 
  draw_plot(boxplots_F1_Levine_32, 0.55, 0.66, 0.4, 0.33) + 
  draw_plot(barplot_mean_F1_pr_re_Levine_32, 0.05, 0.33, 0.4, 0.33) + 
  draw_plot(plot_n_cells_Levine_32, 0.55, 0.36, 0.4, 0.30) + 
  draw_plot(runtime_barplot_Levine_32, 0.05, 0, 0.4, 0.33) + 
  draw_plot(runtime_scatterplot_Levine_32, 0.55, 0.03, 0.4, 0.30) + 
  draw_plot_label(LETTERS[1:6], 
                  c(0, 0.5, 0, 0.5, 0, 0.5), c(0.99, 0.99, 0.66, 0.66, 0.33, 0.33), size = 16)

ggplot2::ggsave("../plots/Levine_2015_marrow_32/plots_multi_panel_Levine2015marrow32.pdf", 
                width = 13, height = 14.5)


# Levine_2015_marrow_13
ggdraw() + 
  draw_plot(barplot_mean_F1_Levine_13, 0.05, 0.66, 0.4, 0.33) + 
  draw_plot(boxplots_F1_Levine_13, 0.55, 0.66, 0.4, 0.33) + 
  draw_plot(barplot_mean_F1_pr_re_Levine_13, 0.05, 0.33, 0.4, 0.33) + 
  draw_plot(plot_n_cells_Levine_13, 0.55, 0.36, 0.4, 0.30) + 
  draw_plot(runtime_barplot_Levine_13, 0.05, 0, 0.4, 0.33) + 
  draw_plot(runtime_scatterplot_Levine_13, 0.55, 0.03, 0.4, 0.30) + 
  draw_plot_label(LETTERS[1:6], 
                  c(0, 0.5, 0, 0.5, 0, 0.5), c(0.99, 0.99, 0.66, 0.66, 0.33, 0.33), size = 16)

ggplot2::ggsave("../plots/Levine_2015_marrow_13/plots_multi_panel_Levine2015marrow13.pdf", 
                width = 13, height = 14.5)


# Nilsson_2013_HSC
ggdraw() + 
  draw_plot(barplot_F1_pr_re_Nilsson, 0.025, 0.64, 0.8, 0.35) + 
  draw_plot(runtime_barplot_Nilsson, 0.025, 0.32, 0.8, 0.32) + 
  draw_plot(runtime_scatterplot_Nilsson, 0.025, 0.03, 0.8, 0.29) + 
  draw_plot_label(LETTERS[1:3], c(0, 0, 0), c(0.99, 0.64, 0.32), size = 16)

ggplot2::ggsave("../plots/Nilsson_2013_HSC/plots_multi_panel_Nilsson2013HSC.pdf", 
                width = 6.5, height = 14.5)


# Mosmann_2014_activ
ggdraw() + 
  draw_plot(barplot_F1_pr_re_Mosmann, 0.025, 0.64, 0.8, 0.35) + 
  draw_plot(runtime_barplot_Mosmann, 0.025, 0.32, 0.8, 0.32) + 
  draw_plot(runtime_scatterplot_Mosmann, 0.025, 0.03, 0.8, 0.29) + 
  draw_plot_label(LETTERS[1:3], c(0, 0, 0), c(0.99, 0.64, 0.32), size = 16)

ggplot2::ggsave("../plots/Mosmann_2014_activ/plots_multi_panel_Mosmann2014activ.pdf", 
                width = 6.5, height = 14.5)




####################
### SESSION INFO ###
####################

# session information for plots

sink(file = "../plots/session_info_plots_main_results.txt")
sessionInfo()
sink()



