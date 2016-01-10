#########################################################################################
# R script to generate additional plots of results (not used in main results sections)
#
# Lukas M. Weber, January 2016
#########################################################################################


### before running this script, run "plots_main_results.R" to load all the necessary packages and objects




###################################################################
### MEAN F1 SCORE: WEIGHTED BY NUMBER OF CELLS PER TRUE CLUSTER ###
###################################################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


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




##########################################################################################
### MEAN F1 SCORE, PRECISION, AND RECALL: WEIGHTED BY NUMBER OF CELLS PER TRUE CLUSTER ###
##########################################################################################

### for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


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

filenames_Levine_32 <- c("../plots/Levine_2015_marrow_32/results_lines_F1_score_by_cluster_size_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_lines_precision_by_cluster_size_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_lines_recall_by_cluster_size_Levine2015marrow32.pdf")
filenames_Levine_13 <- c("../plots/Levine_2015_marrow_13/results_lines_F1_score_by_cluster_size_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_lines_precision_by_cluster_size_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_lines_recall_by_cluster_size_Levine2015marrow13.pdf")


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


