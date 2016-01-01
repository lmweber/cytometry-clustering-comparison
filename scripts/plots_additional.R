#########################################################################################
# R script to generate additional plots of results
#
# Lukas M. Weber, December 2015
#########################################################################################


### before running this script, run the first few sections of "plots_main_results.R" to 
### load all the necessary packages and objects




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


