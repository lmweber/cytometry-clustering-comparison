#########################################################################################
# R script to generate plots of additional results by population
#
# Lukas M. Weber, March 2016
#########################################################################################


### before running this script, run the code in "plots_main_results.R" to load results

library(pheatmap)
library(RColorBrewer)




###########################################################
### HEATMAPS: F1 SCORE, PRECISION, RECALL BY POPULATION ###
###########################################################

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


filenames_Levine_32 <- c("../plots/Levine_2015_marrow_32/by_population/results_by_population_F1_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/by_population/results_by_population_pr_Levine2015marrow32.pdf", 
                         "../plots/Levine_2015_marrow_32/by_population/results_by_population_re_Levine2015marrow32.pdf")

filenames_Levine_13 <- c("../plots/Levine_2015_marrow_13/by_population/results_by_population_F1_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/by_population/results_by_population_pr_Levine2015marrow13.pdf", 
                         "../plots/Levine_2015_marrow_13/by_population/results_by_population_re_Levine2015marrow13.pdf")


# heatmaps of F1 score, precision, and recall for each true population

for (i in 1:length(df_Levine_32)) {
  pheatmap(df_Levine_32[[i]], color = colorRampPalette(brewer.pal(7, "GnBu"))(100), 
           breaks = seq(0, 1, length.out = 100), 
           display_numbers = TRUE, number_color = "black", 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           main = titles_Levine_32[i], filename = filenames_Levine_32[i], 
           width = 6, height = 6)
}

for (i in 1:length(df_Levine_13)) {
  pheatmap(df_Levine_13[[i]], color = colorRampPalette(brewer.pal(7, "GnBu"))(100), 
           breaks = seq(0, 1, length.out = 100), 
           display_numbers = TRUE, number_color = "black", 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           main = titles_Levine_13[i], filename = filenames_Levine_13[i], 
           width = 6, height = 9)
}



