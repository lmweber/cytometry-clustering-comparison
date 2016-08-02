#########################################################################################
# R script to generate main plots of results
#
# Lukas Weber, August 2016
#########################################################################################


library(ggplot2)
library(reshape2)
library(cowplot)  # note masks ggplot2::ggsave()
library(ggrepel)

# helper function for plots
source("../helpers/helper_collapse_data_frame.R")

# load and evaluate results
source("../evaluate_results/evaluate_all_methods.R")
source("../evaluate_results/evaluate_runtime.R")

# color-blind friendly palettes (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cb_pal_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_pal_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ggplot2 default palette (access with ggplot_build(p)$data)
gg_pal <- c("#F8766D", "#00BA38", "#619CFF")




###########################
### PREPARE DATA FRAMES ###
###########################

# results for all methods

res_all <- list(
  #ACCENSE = res_ACCENSE, 
  ClusterX = res_ClusterX, 
  DensVM = res_DensVM, 
  FLOCK = res_FLOCK, 
  flowMeans = res_flowMeans, 
  flowPeaks = res_flowPeaks, 
  FlowSOM = res_FlowSOM, 
  FlowSOM_pre = res_FlowSOM_pre_meta, 
  immunoClust = res_immunoClust, 
  kmeans = res_kmeans, 
  #PhenoGraph = res_PhenoGraph, 
  Rclusterpp = res_Rclusterpp#, 
  #SamSPECTRAL = res_SamSPECTRAL, 
  #SWIFT = res_SWIFT)
)


# collapse into data frames (use helper functions to pad with zeros or NAs for clusters removed by subsampling)

precision_df <- list(
  Levine_32dim = as.data.frame(sapply(res_all, function(res) res[["Levine_32dim"]][["pr"]])), 
  Levine_13dim = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Levine_13dim"]][["pr"]])), 
  Samusik_01   = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Samusik_01"]][["pr"]])), 
  Samusik_all  = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Samusik_all"]][["pr"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["pr"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["pr"]]))
)

precision_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND"]][["pr"]])), 
  FlowCAP_WNV  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV"]][["pr"]])), 
  FlowCAP_ND_alternate  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["pr"]])), 
  FlowCAP_WNV_alternate = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["pr"]]))
)


recall_df <- list(
  Levine_32dim = as.data.frame(sapply(res_all, function(res) res[["Levine_32dim"]][["re"]])), 
  Levine_13dim = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Levine_13dim"]][["re"]])), 
  Samusik_01   = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Samusik_01"]][["re"]])), 
  Samusik_all  = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Samusik_all"]][["re"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["re"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["re"]]))
)

recall_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND"]][["re"]])), 
  FlowCAP_WNV  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV"]][["re"]])), 
  FlowCAP_ND_alternate  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["re"]])), 
  FlowCAP_WNV_alternate = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["re"]]))
)


F1_df <- list(
  Levine_32dim = as.data.frame(sapply(res_all, function(res) res[["Levine_32dim"]][["F1"]])), 
  Levine_13dim = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Levine_13dim"]][["F1"]])), 
  Samusik_01   = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Samusik_01"]][["F1"]])), 
  Samusik_all  = collapse_data_frame_zeros(sapply(res_all, function(res) res[["Samusik_all"]][["F1"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["F1"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["F1"]]))
)

F1_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND"]][["F1"]])), 
  FlowCAP_WNV  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV"]][["F1"]])), 
  FlowCAP_ND_alternate  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["F1"]])), 
  FlowCAP_WNV_alternate = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["F1"]]))
)


labels_matched_df <- list(
  Levine_32dim = as.data.frame(sapply(res_all, function(res) res[["Levine_32dim"]][["labels_matched"]])), 
  Levine_13dim = collapse_data_frame_NAs(sapply(res_all, function(res) res[["Levine_13dim"]][["labels_matched"]])), 
  Samusik_01   = collapse_data_frame_NAs(sapply(res_all, function(res) res[["Samusik_01"]][["labels_matched"]])), 
  Samusik_all  = collapse_data_frame_NAs(sapply(res_all, function(res) res[["Samusik_all"]][["labels_matched"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["labels_matched"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["labels_matched"]]))
)

labels_matched_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND"]][["labels_matched"]])), 
  FlowCAP_WNV  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV"]][["labels_matched"]])), 
  FlowCAP_ND_alternate  = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["labels_matched"]])), 
  FlowCAP_WNV_alternate = as.data.frame(sapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["labels_matched"]]))
)




#####################
### MEAN F1 SCORE ###
#####################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)


# mean F1 score across true populations (unweighted)
mean_F1 <- lapply(F1_df, colMeans)[c("Levine_32dim", "Levine_13dim", "Samusik_01", "Samusik_all")]

# arrange in descending order
mean_F1_ord <- lapply(mean_F1, function(m) {
  ord <- rev(order(m))
  m[ord]
})

# tidy data format (for ggplot)
mean_F1_tidy <- lapply(mean_F1_ord, function(m) {
  d <- data.frame(value = m)
  d["method"] <- factor(rownames(d), levels = rownames(d))
  d
})


# bar plots

for (i in 1:4) {
  nm <- names(mean_F1_tidy)[i]
  title <- paste0("Mean F1 score: ", nm)
  filename <- paste0("../../plots/", nm, "/results_barplot_mean_F1_", nm, ".pdf")
  
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
  
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




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



