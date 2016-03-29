#########################################################################################
# R script to generate plots for ensemble clustering results
#
# Lukas M. Weber, March 2016
#########################################################################################


library(ggplot2)
library(reshape2)
library(cowplot)  # note masks ggplot2::ggsave

# helper function
source("helper_collapse_data_frame.R")

# load ensemble clustering results
source("load_results_ensemble_clustering.R")


### before running this script, also run the code in "plots_main_results.R" to load results




###########################
### PREPARE DATA FRAMES ###
###########################

# combine results into lists

res_Levine_32 <- list(ACCENSE = res_ACCENSE_Levine_32, 
                      DensVM = res_DensVM_Levine_32, 
                      ensemble = res_ensemble_Levine_32, 
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
                      ensemble = res_ensemble_Levine_13, 
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
                    ensemble = res_ensemble_Nilsson, 
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
                    ensemble = res_ensemble_Mosmann, 
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

barplot_mean_F1_Levine_32_ensemble <- 
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

barplot_mean_F1_Levine_32_ensemble
ggplot2::ggsave("../plots/Levine_2015_marrow_32/ensemble_clustering/results_barplot_mean_F1_Levine2015marrow32_ensemble.pdf", 
                width = 5, height = 5)


barplot_mean_F1_Levine_13_ensemble <- 
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

barplot_mean_F1_Levine_13_ensemble
ggplot2::ggsave("../plots/Levine_2015_marrow_13/ensemble_clustering/results_barplot_mean_F1_Levine2015marrow13_ensemble.pdf", 
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

boxplots_F1_Levine_32_ensemble <- 
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

boxplots_F1_Levine_32_ensemble
ggplot2::ggsave("../plots/Levine_2015_marrow_32/ensemble_clustering/results_boxplots_F1_Levine2015marrow32_ensemble.pdf", 
                width = 5, height = 5)


boxplots_F1_Levine_13_ensemble <- 
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

boxplots_F1_Levine_13_ensemble
ggplot2::ggsave("../plots/Levine_2015_marrow_13/ensemble_clustering/results_boxplots_F1_Levine2015marrow13_ensemble.pdf", 
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

barplot_mean_F1_pr_re_Levine_32_ensemble <- 
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

barplot_mean_F1_pr_re_Levine_32_ensemble
ggplot2::ggsave("../plots/Levine_2015_marrow_32/ensemble_clustering/results_barplot_mean_F1_pr_re_Levine2015marrow32_ensemble.pdf", 
                width = 5, height = 5)


barplot_mean_F1_pr_re_Levine_13_ensemble <- 
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

barplot_mean_F1_pr_re_Levine_13_ensemble
ggplot2::ggsave("../plots/Levine_2015_marrow_13/ensemble_clustering/results_barplot_mean_F1_pr_re_Levine2015marrow13_ensemble.pdf", 
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

barplot_F1_pr_re_Nilsson_ensemble <- 
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

barplot_F1_pr_re_Nilsson_ensemble
ggplot2::ggsave("../plots/Nilsson_2013_HSC/ensemble_clustering/results_barplot_F1_pr_re_Nilsson2013HSC_ensemble.pdf", 
                width = 5, height = 5)


barplot_F1_pr_re_Mosmann_ensemble <- 
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

barplot_F1_pr_re_Mosmann_ensemble
ggplot2::ggsave("../plots/Mosmann_2014_activ/ensemble_clustering/results_barplot_F1_pr_re_Mosmann2014activ_ensemble.pdf", 
                width = 5, height = 5)




#########################
### MULTI-PANEL PLOTS ###
#########################

# combine into one multi-panel plot for each data set


# Levine_2015_marrow_32
ggdraw() + 
  draw_plot(barplot_mean_F1_Levine_32_ensemble, 0.05, 0.5, 0.4, 0.5) + 
  draw_plot(boxplots_F1_Levine_32_ensemble, 0.55, 0.5, 0.4, 0.5) + 
  draw_plot(barplot_mean_F1_pr_re_Levine_32_ensemble, 0.05, 0, 0.4, 0.5) + 
  draw_plot_label(LETTERS[1:3], c(0, 0.5, 0), c(1, 1, 0.5), size = 16)

ggplot2::ggsave("../plots/Levine_2015_marrow_32/ensemble_clustering/plots_multi_panel_Levine2015marrow32_ensemble.pdf", 
                width = 13, height = 9.7)


# Levine_2015_marrow_13
ggdraw() + 
  draw_plot(barplot_mean_F1_Levine_13_ensemble, 0.05, 0.5, 0.4, 0.5) + 
  draw_plot(boxplots_F1_Levine_13_ensemble, 0.55, 0.5, 0.4, 0.5) + 
  draw_plot(barplot_mean_F1_pr_re_Levine_13_ensemble, 0.05, 0, 0.4, 0.5) + 
  draw_plot_label(LETTERS[1:3], c(0, 0.5, 0), c(1, 1, 0.5), size = 16)

ggplot2::ggsave("../plots/Levine_2015_marrow_13/ensemble_clustering/plots_multi_panel_Levine2015marrow13_ensemble.pdf", 
                width = 13, height = 9.7)


# Nilsson_2013_HSC
ggdraw() + 
  draw_plot(barplot_F1_pr_re_Nilsson_ensemble, 0.3, 0, 0.4, 1)

ggplot2::ggsave("../plots/Nilsson_2013_HSC/ensemble_clustering/plots_multi_panel_Nilsson2013HSC_ensemble.pdf", 
                width = 13, height = 4.9)


# Mosmann_2014_activ
ggdraw() + 
  draw_plot(barplot_F1_pr_re_Mosmann_ensemble, 0.3, 0, 0.4, 1)

ggplot2::ggsave("../plots/Mosmann_2014_activ/ensemble_clustering/plots_multi_panel_Mosmann2014activ_ensemble.pdf", 
                width = 13, height = 4.9)



