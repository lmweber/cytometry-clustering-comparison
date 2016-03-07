#########################################################################################
# R script to generate plots of results weighted by population size
#
# Lukas M. Weber, March 2016
#########################################################################################


### before running this script, run the code in "plots_main_results.R" to load results



##################################################
### MEAN F1 SCORE: WEIGHTED BY POPULATION SIZE ###
##################################################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


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

barplot_mean_F1_Levine_32_weighted <- 
  ggplot(mean_F1_Levine_32_tidy_weighted, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.3f", round(value, 3)), y = value + 0.07, angle = 90), size = 3.5) + 
  ylim(0, 1.02) + 
  ylab("mean F1 score (weighted)") + 
  ggtitle("Mean F1 score (weighted): Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

barplot_mean_F1_Levine_32_weighted
ggplot2::ggsave("../plots/Levine_2015_marrow_32/weighted_means/results_mean_F1_weighted_Levine2015marrow32.pdf", 
                height = 5, width = 5)


barplot_mean_F1_Levine_13_weighted <- 
  ggplot(mean_F1_Levine_13_tidy_weighted, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.3f", round(value, 3)), y = value + 0.07, angle = 90), size = 3.5) + 
  ylim(0, 1.02) + 
  ylab("mean F1 score (weighted)") + 
  ggtitle("Mean F1 score (weighted): Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

barplot_mean_F1_Levine_13_weighted
ggplot2::ggsave("../plots/Levine_2015_marrow_13/weighted_means/results_mean_F1_weighted_Levine2015marrow13.pdf", 
                height = 5, width = 5)




#####################################################################
### MEAN F1 SCORE, PRECISION, RECALL: WEIGHTED BY POPULATION SIZE ###
#####################################################################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


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

barplot_mean_F1_pr_re_Levine_32_weighted <- 
  ggplot(plot_data_Levine_32_weighted, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1.02) + 
  ylab("") + 
  ggtitle("Mean F1, precision, recall (weighted): Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 11), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.75, 0.96), 
        legend.direction = "horizontal", 
        legend.key.size = unit(3, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

barplot_mean_F1_pr_re_Levine_32_weighted
ggplot2::ggsave("../plots/Levine_2015_marrow_32/weighted_means/results_mean_F1_pr_re_weighted_Levine2015marrow32.pdf", 
                height = 5, width = 5)


barplot_mean_F1_pr_re_Levine_13_weighted <- 
  ggplot(plot_data_Levine_13_weighted, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1.02) + 
  ylab("") + 
  ggtitle("Mean F1, precision, recall (weighted): Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 11), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.75, 0.96), 
        legend.direction = "horizontal", 
        legend.key.size = unit(3, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

barplot_mean_F1_pr_re_Levine_13_weighted
ggplot2::ggsave("../plots/Levine_2015_marrow_13/weighted_means/results_mean_F1_pr_re_weighted_Levine2015marrow13.pdf", 
                height = 5, width = 5)




#########################
### MULTI-PANEL PLOTS ###
#########################

# combine into one multi-panel plot for each data set


# Levine_32
ggdraw() + 
  draw_plot(barplot_mean_F1_Levine_32_weighted, 0.05, 0, 0.4, 1) + 
  draw_plot(barplot_mean_F1_pr_re_Levine_32_weighted, 0.55, 0, 0.4, 1) +
  draw_plot_label(LETTERS[1:2], c(0, 0.5), c(1, 1), size = 16)

ggplot2::ggsave("../plots/Levine_2015_marrow_32/weighted_means/weighted_means_multi_panel_Levine2015marrow32.pdf", 
                width = 13, height = 5)


# Levine_13
ggdraw() + 
  draw_plot(barplot_mean_F1_Levine_13_weighted, 0.05, 0, 0.4, 1) + 
  draw_plot(barplot_mean_F1_pr_re_Levine_13_weighted, 0.55, 0, 0.4, 1) +
  draw_plot_label(LETTERS[1:2], c(0, 0.5), c(1, 1), size = 16)

ggplot2::ggsave("../plots/Levine_2015_marrow_13/weighted_means/weighted_means_multi_panel_Levine2015marrow13.pdf", 
                width = 13, height = 5)



