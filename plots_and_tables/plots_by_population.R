#########################################################################################
# R script to generate plots of additional results by population
#
# Lukas Weber, September 2016
#########################################################################################


library(pheatmap)
library(RColorBrewer)

# load main results (from main plots file)
load("main_results.RData")




###########################################################
### HEATMAPS: F1 SCORE, PRECISION, RECALL BY POPULATION ###
###########################################################

# only for data sets with multiple populations (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)


# rows (true populations) arranged by size; columns (methods) arranged in the same order as previously

ord_size <- lapply(tbl_truth_multiple, function(tb) rev(order(tb)))

plot_data_by_pop <- vector("list", length(ord_size))
names(plot_data_by_pop) <- names(ord_size)

for (i in 1:length(plot_data_by_pop)) {
  plot_data_by_pop[[i]] <- list(F1 = F1_df[[i]][ord_size[[i]], ord[[i]]], 
                                pr = precision_df[[i]][ord_size[[i]], ord[[i]]], 
                                re = recall_df[[i]][ord_size[[i]], ord[[i]]])
}


# heatmaps

title_prefixes <- paste0("(", LETTERS[1:3], ")   ", c("F1 score: ", "Precision: ", "Recall: "))
fig_height <- c(6, 9, 9, 9)

for (i in 1:length(plot_data_by_pop)) {  ## loop over data sets
  for (j in 1:length(plot_data_by_pop[[i]])) {  ## loop over F1, pr, re
    
    nm_dataset <- names(plot_data_by_pop)[i]
    nm_F1_pr_re <- names(plot_data_by_pop[[i]])[j]
    
    pheatmap(plot_data_by_pop[[i]][[j]], 
             color = colorRampPalette(brewer.pal(7, "GnBu"))(100), 
             breaks = seq(0, 1, length.out = 100), 
             display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
             cluster_rows = FALSE, cluster_cols = FALSE, 
             main = paste0(title_prefixes[j], nm_dataset), 
             filename = paste0("../../plots/", nm_dataset, "/by_population/heatmap_by_population_", 
                               nm_F1_pr_re, "_", nm_dataset, ".pdf"), 
             width = 6, height = fig_height[i])
  }
}



