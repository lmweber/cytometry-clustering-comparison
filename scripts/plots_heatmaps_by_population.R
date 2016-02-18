#########################################################################################
# R script to generate heatmaps with additional results by cell population
#
# Lukas M. Weber, February 2016
#########################################################################################


### before running this script, run the code in "plots_main_results.R" to load results

library(pheatmap)
library(RColorBrewer)



#################################################################################
### HEATMAPS: F1 SCORE, PRECISION, AND RECALL BY POPULATION (ORDERED BY SIZE) ###
#################################################################################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)


# rows (true populations) ranked by size, and columns (methods) arranged in the same order as previously

ord_size_rev_Levine_32 <- rev(order(tbl_truth_Levine_32))
ord_size_rev_Levine_13 <- rev(order(tbl_truth_Levine_13))


df_Levine_32 <- list(F1_df_Levine_32[ord_size_rev_Levine_32, ord_Levine_32], 
                     precision_df_Levine_32[ord_size_rev_Levine_32, ord_Levine_32], 
                     recall_df_Levine_32[ord_size_rev_Levine_32, ord_Levine_32])

df_Levine_13 <- list(F1_df_Levine_13[ord_size_rev_Levine_13, ord_Levine_13], 
                     precision_df_Levine_13[ord_size_rev_Levine_13, ord_Levine_13], 
                     recall_df_Levine_13[ord_size_rev_Levine_13, ord_Levine_13])


titles_Levine_32 <- c("(A)   F1 score: Levine_2015_marrow_32", 
                      "(B)   Precision: Levine_2015_marrow_32", 
                      "(C)   Recall: Levine_2015_marrow_32")

titles_Levine_13 <- c("(A)   F1 score: Levine_2015_marrow_13", 
                      "(B)   Precision: Levine_2015_marrow_13", 
                      "(C)   Recall: Levine_2015_marrow_13")


filenames_Levine_32 <- c("../plots/Levine_2015_marrow_32/results_heatmap_F1_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_heatmap_precision_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/results_heatmap_recall_Levine2015marrow32.pdf")

filenames_Levine_13 <- c("../plots/Levine_2015_marrow_13/results_heatmap_F1_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_heatmap_precision_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/results_heatmap_recall_Levine2015marrow13.pdf")


# heatmaps of F1 score, precision, and recall for each true population

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



