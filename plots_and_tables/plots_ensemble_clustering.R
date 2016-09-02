#########################################################################################
# R script to generate plots for ensemble clustering results
#
# Lukas Weber, September 2016
#########################################################################################


library(ggplot2)
library(reshape2)
library(cowplot)  # note masks ggplot2::ggsave()

# load main results (from main plots file)
load("main_results.RData")

# load ensemble clustering results
source("../ensemble_clustering/evaluate_ensemble_clustering.R")




###########################
### PREPARE DATA FRAMES ###
###########################

# results for all methods

res_all <- list(
  ACCENSE = res_ACCENSE, 
  ClusterX = res_ClusterX, 
  DensVM = res_DensVM, 
  ensemble = res_ensemble, 
  FLOCK = res_FLOCK, 
  flowClust = res_flowClust, 
  flowMeans = res_flowMeans, 
  flowMerge = res_flowMerge, 
  flowPeaks = res_flowPeaks, 
  FlowSOM = res_FlowSOM, 
  FlowSOM_pre = res_FlowSOM_pre, 
  immunoClust = res_immunoClust, 
  kmeans = res_kmeans, 
  PhenoGraph = res_PhenoGraph, 
  Rclusterpp = res_Rclusterpp, 
  SamSPECTRAL = res_SamSPECTRAL, 
  SPADE = res_SPADE, 
  SWIFT = res_SWIFT, 
  Xshift = res_Xshift
)


# collapse into data frames (use helper functions to pad with zeros or NAs for missing
# methods or true populations removed by subsampling)

precision_df <- list(
  Levine_32dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_32dim"]][["pr"]])), 
  Levine_13dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_13dim"]][["pr"]])), 
  Samusik_01   = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_01"]][["pr"]])), 
  Samusik_all  = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_all"]][["pr"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["pr"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["pr"]]))
)

recall_df <- list(
  Levine_32dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_32dim"]][["re"]])), 
  Levine_13dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_13dim"]][["re"]])), 
  Samusik_01   = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_01"]][["re"]])), 
  Samusik_all  = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_all"]][["re"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["re"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["re"]]))
)

F1_df <- list(
  Levine_32dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_32dim"]][["F1"]])), 
  Levine_13dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_13dim"]][["F1"]])), 
  Samusik_01   = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_01"]][["F1"]])), 
  Samusik_all  = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_all"]][["F1"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["F1"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["F1"]]))
)




############################
### PLOTS: MEAN F1 SCORE ###
############################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

# mean F1 score across true populations (unweighted)
mean_F1 <- lapply(F1_df, colMeans)[data_sets_multiple]

# arrange in descending order
ord <- lapply(mean_F1, function(m) rev(order(m)))
for (i in 1:length(mean_F1)) {
  mean_F1[[i]] <- mean_F1[[i]][ord[[i]]]
}

# tidy data format (for ggplot)
mean_F1_tidy <- lapply(mean_F1, function(m) {
  d <- data.frame(value = m)
  d["method"] <- factor(rownames(d), levels = rownames(d))
  d
})

# bar plots

barplots_mean_F1 <- vector("list", length(mean_F1_tidy))
names(barplots_mean_F1) <- names(mean_F1_tidy)

for (i in 1:4) {
  nm <- names(mean_F1_tidy)[i]
  title <- paste0("Mean F1 score: ", nm)
  filename <- paste0("../../plots/", nm, "/ensemble_clustering/results_barplot_mean_F1_ensemble_", nm, ".pdf")
  
  pl <- 
    ggplot(mean_F1_tidy[[i]], aes(x = method, y = value)) + 
    geom_bar(stat = "identity", fill = "royalblue3") + 
    geom_text(aes(label = sprintf("%.3f", round(value, 3)), y = value + 0.08, angle = 90), size = 3.5) + 
    ylim(0, 1) + 
    ylab("mean F1 score") + 
    ggtitle(title) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  barplots_mean_F1[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




#################################
### PLOTS: F1 SCORE BOX PLOTS ###
#################################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

# arrange in same order as previous plots
F1_df_multiple <- F1_df[data_sets_multiple]
for (i in 1:length(F1_df_multiple)) {
  F1_df_multiple[[i]] <- F1_df_multiple[[i]][, ord[[i]]]
}

# tidy data format (for ggplot)
F1_df_tidy <- lapply(F1_df_multiple, function(m) {
  d <- data.frame(value = as.vector(as.matrix(m)))
  d["method"] <- rep(factor(colnames(m), levels = colnames(m)), each = nrow(m))
  d
})

# box plots

boxplots_F1 <- vector("list", length(F1_df_tidy))
names(boxplots_F1) <- names(F1_df_tidy)

for (i in 1:4) {
  nm <- names(F1_df_tidy)[i]
  title <- paste0("F1 score: ", nm)
  filename <- paste0("../../plots/", nm, "/ensemble_clustering/results_boxplots_F1_ensemble_", nm, ".pdf")
  
  pl <- 
    ggplot(F1_df_tidy[[i]], aes(x = method, y = value)) + 
    geom_boxplot(col = "gray50", fill = "aliceblue") + 
    geom_point(shape = 1, col = "darkblue") + 
    stat_summary(fun.y = mean, color = "red", geom = "point", shape = 1) + 
    ylim(0, 1) + 
    ylab("F1 score") + 
    ggtitle(title) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  boxplots_F1[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




#########################################################
### PLOTS: MEAN F1 SCORE, MEAN PRECISION, MEAN RECALL ###
#########################################################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

# mean precision and mean recall across true populations (unweighted)
mean_precision <- lapply(precision_df, colMeans)[data_sets_multiple]
mean_recall <- lapply(recall_df, colMeans)[data_sets_multiple]

# arrange in same order as previous plots
for (i in 1:length(mean_precision)) {
  mean_precision[[i]] <- mean_precision[[i]][ord[[i]]]
  mean_recall[[i]] <- mean_recall[[i]][ord[[i]]]
}

# tidy data format (for ggplot)
f_plot_data <- function(f, p, r) {
  d <- data.frame(F1_score = f, precision = p, recall = r)
  d["method"] <- factor(rownames(d), levels = rownames(d))
  d <- melt(d, id.vars = "method", measure.vars = c("F1_score", "precision", "recall"))
  d
}
plot_data <- mapply(f_plot_data, mean_F1, mean_precision, mean_recall, SIMPLIFY = FALSE)

# bar plots of mean F1 score, mean precision, mean recall (in same order as previously)

barplots_mean_F1_pr_re <- vector("list", length(plot_data))
names(barplots_mean_F1_pr_re) <- names(plot_data)

for (i in 1:4) {
  nm <- names(plot_data)[i]
  title <- paste0("Mean F1 score, precision, recall: ", nm)
  filename <- paste0("../../plots/", nm, "/ensemble_clustering/results_barplot_mean_F1_pr_re_ensemble_", nm, ".pdf")
  
  pl <- 
    ggplot(plot_data[[i]], aes(x = method, y = value, group = variable, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
    ylim(0, 1) + 
    ylab("") + 
    ggtitle(title) + 
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
  
  barplots_mean_F1_pr_re[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




############################################################
### PLOTS: RARE POPULATIONS: F1 SCORE, PRECISION, RECALL ###
############################################################

# for data sets with a single rare population of interest (Nilsson_rare, Mosmann_rare)

# arrange by decreasing F1 score
ord_rare <- lapply(F1_df[data_sets_single], function(d) rev(order(d)))

precision_rare <- precision_df[data_sets_single]
recall_rare <- recall_df[data_sets_single]
F1_rare <- F1_df[data_sets_single]

for (i in 1:length(precision_rare)) {
  precision_rare[[i]] <- unlist(precision_rare[[i]][ord_rare[[i]]])
  recall_rare[[i]] <- unlist(recall_rare[[i]][ord_rare[[i]]])
  F1_rare[[i]] <- unlist(F1_rare[[i]][ord_rare[[i]]])
}

# tidy data format (for ggplot)
plot_data_rare <- mapply(f_plot_data, F1_rare, precision_rare, recall_rare, SIMPLIFY = FALSE)

# bar plots

barplots_F1_pr_re <- vector("list", length(plot_data_rare))
names(barplots_F1_pr_re) <- names(plot_data_rare)

for (i in 1:2) {
  nm <- names(plot_data_rare)[i]
  title <- paste0("Rare population: ", nm)
  filename <- paste0("../../plots/", nm, "/ensemble_clustering/results_barplot_F1_pr_re_ensemble_", nm, ".pdf")
  
  pl <- 
    ggplot(plot_data_rare[[i]], aes(x = method, y = value, group = variable, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
    ylim(0, 1.05) + 
    ylab("") + 
    ggtitle(title) + 
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
  
  barplots_F1_pr_re[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




##########################
### PLOTS: MULTI-PANEL ###
##########################

# combine into one multi-panel plot for each data set

# data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

multi_panel_multiple <- vector("list", length(mean_F1_tidy))
names(multi_panel_multiple) <- names(mean_F1_tidy)

for (i in 1:4) {
  nm <- names(multi_panel_multiple)[i]
  filename <- paste0("../../plots/", nm, "/ensemble_clustering/plots_multi_panel_ensemble_", nm, ".pdf")
  
  pl <- ggdraw() + 
    draw_plot(barplots_mean_F1[[i]], 0.05, 0.5, 0.4, 0.5) + 
    draw_plot(boxplots_F1[[i]], 0.55, 0.5, 0.4, 0.5) + 
    draw_plot(barplots_mean_F1_pr_re[[i]], 0.05, 0, 0.4, 0.5) + 
    draw_plot_label(LETTERS[1:3], 
                    c(0, 0.5, 0), c(0.99, 0.99, 0.5), size = 16)
  
  multi_panel_multiple[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, width = 13, height = 9.7)
}


## multi-panel plots not required for data sets with a single rare population of interest
## (Nilsson_rare, Mosmann_rare), since there is only one panel



