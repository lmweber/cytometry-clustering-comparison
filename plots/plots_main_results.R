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
CURRENT_DIR <- getwd()
setwd("../evaluate_results")
source("evaluate_all_methods.R")  ## takes 15 min
source("evaluate_runtime.R")
setwd(CURRENT_DIR)

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
  FlowSOM_pre = res_FlowSOM_pre, 
  immunoClust = res_immunoClust, 
  kmeans = res_kmeans, 
  #PhenoGraph = res_PhenoGraph, 
  Rclusterpp = res_Rclusterpp, 
  SamSPECTRAL = res_SamSPECTRAL#, 
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


# population sizes: number of assigned cells per true population

files_truth <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik/data/Samusik_all.fcs"), 
  Nilsson_rare = file.path(DATA_DIR, "Nilsson_rare/data/Nilsson_rare.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare.fcs"), 
  FlowCAP_ND   = file.path(DATA_DIR, "FlowCAP_ND/data/FlowCAP_ND.fcs"), 
  FlowCAP_WNV  = file.path(DATA_DIR, "FlowCAP_WNV/data/FlowCAP_WNV.fcs")
)

clus_truth <- lapply(files_truth, function(f) {
  d <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
  d[, "label"]
})

sapply(clus_truth, length)
tbl_truth <- lapply(clus_truth, table)
tbl_truth
sapply(tbl_truth, length)




#########################
### TABLES OF RESULTS ###
#########################

data_sets_multiple <- 1:4
data_sets_single   <- 5:6
data_sets_FlowCAP  <- 7:8
data_sets_FlowCAP_alternate <- 9:10


# mean F1 score (multiple populations) or F1 score (single population) for each method and data set

# data sets with multiple populations (mean F1 score)
sapply(res_all, function(r) sapply(r[data_sets_multiple], function(s) s$mean_F1))
# data sets with single population of interest (F1 score)
sapply(res_all, function(r) sapply(r[data_sets_single], function(s) s$F1))


# FlowCAP data sets (Hungarian algorithm matching, unweighted averages)
sapply(res_all, function(r) sapply(r[data_sets_FlowCAP], function(s) s$mean_F1))

# FlowCAP data sets: alternate (max F1 score matching, averages weighted by number of cells)
sapply(res_all, function(r) sapply(r[data_sets_FlowCAP_alternate], function(s) s$mean_F1))




#####################
### MEAN F1 SCORE ###
#####################

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

for (i in 1:4) {
  nm <- names(F1_df_tidy)[i]
  title <- paste0("F1 score: ", nm)
  filename <- paste0("../../plots/", nm, "/results_boxplots_F1_", nm, ".pdf")
  
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
  
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




##################################################
### MEAN F1 SCORE, MEAN PRECISION, MEAN RECALL ###
##################################################

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

for (i in 1:4) {
  nm <- names(plot_data)[i]
  title <- paste0("Mean F1 score, precision, recall: ", nm)
  filename <- paste0("../../plots/", nm, "/results_barplot_mean_F1_pr_re_", nm, ".pdf")
  
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
  
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




########################
### POPULATION SIZES ###
########################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

# tidy data format (for ggplot)
tbl_truth_multiple <- tbl_truth[data_sets_multiple]
n_cells_tidy <- lapply(tbl_truth_multiple, function(d) {
  d <- as.data.frame(d)
  colnames(d) <- c("population", "value")
  d
})

# plots of population sizes: number of assigned cells per true population

ymaxs <- list(30000, 15750, 15500, 132500)
offsets <- list(600, 300, 300, 2500)

for (i in 1:4) {
  nm <- names(n_cells_tidy)[i]
  title <- paste0("Manually gated populations: ", nm)
  filename <- paste0("../../plots/", nm, "/results_no_of_cells_", nm, ".pdf")
  
  pl <- 
    ggplot(n_cells_tidy[[i]], aes(x = population, y = value)) + 
    geom_bar(stat = "identity", fill = "darkgray") + 
    geom_text(aes(label = value, y = value + offsets[[i]], angle = 90), hjust = "left", size = 3.5) + 
    scale_y_continuous(limits = c(0, ymaxs[[i]]), labels = scales::comma) + 
    xlab("manually gated population") + 
    ylab("number of cells") + 
    ggtitle(title) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.text.x = element_text(size = 9))
  
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




#####################################################
### RARE POPULATIONS: F1 SCORE, PRECISION, RECALL ###
#####################################################

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

for (i in 1:2) {
  nm <- names(plot_data_rare)[i]
  title <- paste0("Rare population: ", nm)
  filename <- paste0("../../plots/", nm, "/results_barplot_F1_pr_re_", nm, ".pdf")
  
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
  
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




##########################
### RUNTIME: BAR PLOTS ###
##########################

# all data sets

# arrange runtimes in ascending order
res_runtime_ord <- lapply(res_runtime, function(r) {
  r_ord <- unlist(r)
  r_ord[order(r_ord)]
})

# tidy data format (for ggplot)
runtime_tidy <- lapply(res_runtime_ord, function(r) {
  d <- data.frame(value = r)
  d["method"] <- factor(rownames(d), levels = rownames(d))
  d
})

# which methods required subsampling
runtime_tidy[["Levine_32dim"]]["subsampling"] <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*")
runtime_tidy[["Levine_13dim"]]["subsampling"] <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*")
runtime_tidy[["Samusik_01"]]["subsampling"]   <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*")
runtime_tidy[["Samusik_all"]]["subsampling"]  <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*")
runtime_tidy[["Nilsson_rare"]]["subsampling"] <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*")
runtime_tidy[["Mosmann_rare"]]["subsampling"] <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*")

# which methods required multiple cores
runtime_tidy[["Levine_32dim"]]["cores"] <- c("^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^")
runtime_tidy[["Levine_13dim"]]["cores"] <- c("^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^")
runtime_tidy[["Samusik_01"]]["cores"]   <- c("^", "^", "^", "^", "^", "^", "^", "^", "^", "^")
runtime_tidy[["Samusik_all"]]["cores"]  <- c("^", "^", "^", "^", "^", "^", "^", "^", "^", "^")
runtime_tidy[["Nilsson_rare"]]["cores"] <- c("^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^")
runtime_tidy[["Mosmann_rare"]]["cores"] <- c("^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^", "^")


# bar plots

offsets <- list(1000, 200, 500, 600, 150, 600)
ymaxs <- list(50000, 11000, 25000, 29000, 6500, 30000)

for (i in 1:6) {
  nm <- names(runtime_tidy)[i]
  title <- paste0("Runtime: ", nm)
  filename <- paste0("../../plots/", nm, "/runtime_barplot_", nm, ".pdf")
  
  pl <- 
    ggplot(runtime_tidy[[i]], aes(x = method, y = value)) + 
    geom_bar(stat = "identity", fill = "mediumpurple") + 
    geom_text(aes(label = paste0(round(value, 0), " ", subsampling, cores), 
                  y = value + offsets[[i]], angle = 90), hjust = "left", size = 3.5) + 
    ggtitle(title) + 
    ylim(0, ymaxs[[i]]) + 
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
  
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




##########################################################
### RUNTIME: SCATTER PLOTS (RUNTIME VS. MEAN F1 SCORE) ###
##########################################################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)


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



