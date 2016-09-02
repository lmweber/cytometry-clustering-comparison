#########################################################################################
# R script to generate plots and tables of main results
#
# Lukas Weber, September 2016
#########################################################################################


library(ggplot2)
library(reshape2)
library(cowplot)  # note masks ggplot2::ggsave()
library(ggrepel)
library(lubridate)
library(xtable)

# helper function for plots
source("../helpers/helper_collapse_df.R")

# load and evaluate results
CURRENT_DIR <- getwd()
setwd("../evaluate_results")
source("evaluate_all_methods.R")  ## runtime: 20 min
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
  ACCENSE = res_ACCENSE, 
  ClusterX = res_ClusterX, 
  DensVM = res_DensVM, 
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

precision_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_ND"]][["mean_pr"]])))), 
  FlowCAP_WNV  = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_WNV"]][["mean_pr"]])))), 
  FlowCAP_ND_alternate  = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["mean_pr"]])))), 
  FlowCAP_WNV_alternate = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["mean_pr"]]))))
)


recall_df <- list(
  Levine_32dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_32dim"]][["re"]])), 
  Levine_13dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_13dim"]][["re"]])), 
  Samusik_01   = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_01"]][["re"]])), 
  Samusik_all  = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_all"]][["re"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["re"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["re"]]))
)

recall_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_ND"]][["mean_re"]])))), 
  FlowCAP_WNV  = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_WNV"]][["mean_re"]])))), 
  FlowCAP_ND_alternate  = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["mean_re"]])))), 
  FlowCAP_WNV_alternate = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["mean_re"]]))))
)


F1_df <- list(
  Levine_32dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_32dim"]][["F1"]])), 
  Levine_13dim = collapse_df_zeros(lapply(res_all, function(res) res[["Levine_13dim"]][["F1"]])), 
  Samusik_01   = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_01"]][["F1"]])), 
  Samusik_all  = collapse_df_zeros(lapply(res_all, function(res) res[["Samusik_all"]][["F1"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["F1"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["F1"]]))
)

F1_df_FlowCAP <- list(
  FlowCAP_ND   = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_ND"]][["mean_F1"]])))), 
  FlowCAP_WNV  = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_WNV"]][["mean_F1"]])))), 
  FlowCAP_ND_alternate  = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_ND_alternate"]][["mean_F1"]])))), 
  FlowCAP_WNV_alternate = as.data.frame(t(unlist(lapply(res_all, function(res) res[["FlowCAP_WNV_alternate"]][["mean_F1"]]))))
)


labels_matched_df <- list(
  Levine_32dim = collapse_df_NAs(lapply(res_all, function(res) res[["Levine_32dim"]][["labels_matched"]])), 
  Levine_13dim = collapse_df_NAs(lapply(res_all, function(res) res[["Levine_13dim"]][["labels_matched"]])), 
  Samusik_01   = collapse_df_NAs(lapply(res_all, function(res) res[["Samusik_01"]][["labels_matched"]])), 
  Samusik_all  = collapse_df_NAs(lapply(res_all, function(res) res[["Samusik_all"]][["labels_matched"]])), 
  Nilsson_rare = as.data.frame(lapply(res_all, function(res) res[["Nilsson_rare"]][["labels_matched"]])), 
  Mosmann_rare = as.data.frame(lapply(res_all, function(res) res[["Mosmann_rare"]][["labels_matched"]]))
)

# note: labels_matched not available for FlowCAP data sets



# ---------------------
# true population sizes
# ---------------------

# number of assigned (i.e. manually gated) cells per true population

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

n_cells_truth <- tbl_truth
n_cells_truth

sapply(tbl_truth, length)  # number of true populations




##############################
### MAIN TABLES OF RESULTS ###
##############################

data_sets_multiple <- 1:4
data_sets_single   <- 5:6
data_sets_FlowCAP  <- 7:8
data_sets_FlowCAP_alternate <- 9:10


# --------------------------
# mean F1 scores / F1 scores
# --------------------------

# data sets with multiple populations of interest (mean F1 score)
tbl_multiple <- sapply(res_all, function(r) {
  sapply(r[data_sets_multiple], function(s) {
    ifelse(is.null(s$mean_F1), NA, round(as.numeric(s$mean_F1), 3))
  })
})
tbl_multiple

# data sets with single rare population of interest (F1 score)
tbl_single <- sapply(res_all, function(r) {
  sapply(r[data_sets_single], function(s) {
    ifelse(is.null(s$F1), NA, round(as.numeric(s$F1), 3))
  })
})
tbl_single

# combined table with formatting
tbl_combined <- as.data.frame(cbind(t(tbl_multiple), t(tbl_single)))
tbl_combined



# --------
# runtimes
# --------

tbl_runtime <- as.data.frame(sapply(res_runtime[1:6], as.data.frame))

# re-order
tbl_runtime <- tbl_runtime[rownames(tbl_combined), ]
tbl_runtime

tbl_runtime_formatted <- sapply(tbl_runtime, function(col) {
  time <- round(as.numeric(col), 0)
  time <- lubridate::seconds_to_period(time)
  sprintf("%02d:%02d:%02d", lubridate::hour(time), lubridate::minute(time), lubridate::second(time))
})
rownames(tbl_runtime_formatted) <- rownames(tbl_runtime)
colnames(tbl_runtime_formatted) <- paste0(colnames(tbl_runtime_formatted), "_time")
tbl_runtime_formatted <- as.data.frame(tbl_runtime_formatted)
tbl_runtime_formatted



# --------------------------
# combined (with formatting)
# --------------------------

tbl_main <- cbind(tbl_combined, tbl_runtime_formatted)[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]
tbl_main

# export table in Latex format
sink("../../tables/table_main_results.txt")
print(xtable(tbl_main, digits = 3))
sink()



# -----------------
# FlowCAP data sets
# -----------------

# FlowCAP data sets (Hungarian algorithm for cluster matching; unweighted averages)
tbl_FlowCAP <- sapply(res_all, function(r) {
  sapply(r[data_sets_FlowCAP], function(s) {
    ifelse(is.null(s$mean_F1), NA, round(as.numeric(s$mean_F1), 3))
  })
})
rownames(tbl_FlowCAP) <- c("FlowCAP_ND", "FlowCAP_WNV")
tbl_FlowCAP <- as.data.frame(t(tbl_FlowCAP))

tbl_FlowCAP


# FlowCAP data sets: alternate (maximum F1 score for cluster matching; averages weighted by number of cells)
tbl_FlowCAP_alt <- sapply(res_all, function(r) {
  sapply(r[data_sets_FlowCAP_alternate], function(s) {
    ifelse(is.null(s$mean_F1), NA, round(as.numeric(s$mean_F1), 3))
  })
})
rownames(tbl_FlowCAP_alt) <- c("FlowCAP_ND_alt", "FlowCAP_WNV_alt")
tbl_FlowCAP_alt <- as.data.frame(t(tbl_FlowCAP_alt))

tbl_FlowCAP_alt


# combined table (with formatting); export in Latex format
tbl_FlowCAP_combined <- cbind(tbl_FlowCAP, " " = "", tbl_FlowCAP_alt)
tbl_FlowCAP_combined

# export table in Latex format
sink("../../tables/table_FlowCAP_results.txt")
print(xtable(tbl_FlowCAP_combined, digits = 3))
sink()



# ---------------------------------------
# alternative calculations (for checking)
# ---------------------------------------

# data sets with multiple populations
sapply(res_all, function(r) sapply(r[data_sets_multiple], function(s) s$mean_F1))
lapply(lapply(F1_df[data_sets_multiple], colMeans), t)

# data sets with single population
sapply(res_all, function(r) sapply(r[data_sets_single], function(s) s$F1))
F1_df[data_sets_single]

# FlowCAP data sets
sapply(res_all, function(r) sapply(r[data_sets_FlowCAP], function(s) s$mean_F1))
F1_df_FlowCAP[1:2]

# FlowCAP data sets: alternate
sapply(res_all, function(r) sapply(r[data_sets_FlowCAP_alternate], function(s) s$mean_F1))
F1_df_FlowCAP[3:4]




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
  
  barplots_mean_F1_pr_re[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




###############################
### PLOTS: POPULATION SIZES ###
###############################

# for data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

# tidy data format (for ggplot)
tbl_truth_multiple <- tbl_truth[data_sets_multiple]
n_cells_tidy <- lapply(tbl_truth_multiple, function(d) {
  d <- as.data.frame(d)
  colnames(d) <- c("population", "value")
  d
})

# plots of population sizes: number of assigned cells per true population

no_of_cells <- vector("list", length(n_cells_tidy))
names(no_of_cells) <- names(n_cells_tidy)

ymaxs_ncells <- list(30000, 15750, 15500, 132500)
offsets_ncells <- list(600, 300, 300, 2500)

for (i in 1:4) {
  nm <- names(n_cells_tidy)[i]
  title <- paste0("Manually gated populations: ", nm)
  filename <- paste0("../../plots/", nm, "/results_no_of_cells_", nm, ".pdf")
  
  pl <- 
    ggplot(n_cells_tidy[[i]], aes(x = population, y = value)) + 
    geom_bar(stat = "identity", fill = "darkgray") + 
    geom_text(aes(label = value, y = value + offsets_ncells[[i]], angle = 90), hjust = "left", size = 3.5) + 
    scale_y_continuous(limits = c(0, ymaxs_ncells[[i]]), labels = scales::comma) + 
    xlab("manually gated population") + 
    ylab("number of cells") + 
    ggtitle(title) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.text.x = element_text(size = 9))
  
  no_of_cells[[i]] <- pl
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
  
  barplots_F1_pr_re[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




################################
### PLOTS: RUNTIME BAR PLOTS ###
################################

# arrange runtimes in ascending order
res_runtime_ord <- lapply(res_runtime, function(r) {
  r <- r[!is.na(r)]
  r_ord <- unlist(r)
  r_ord[order(r_ord)]
})

# tidy data format (for ggplot)
runtime_tidy <- lapply(res_runtime_ord[c(data_sets_multiple, data_sets_single)], function(r) {
  d <- data.frame(value = r)
  d["method"] <- factor(rownames(d), levels = rownames(d))
  d
})

# plotting symbols
symbol_sub <- "*"
symbol_cores <- "^"

# which methods required subsampling (from script "evaluate_runtime.R")
runtime_tidy[["Levine_32dim"]][, "subsampling"] <- ""
runtime_tidy[["Levine_32dim"]][which_sub_Levine_32dim, "subsampling"] <- symbol_sub
runtime_tidy[["Levine_13dim"]][, "subsampling"] <- ""
runtime_tidy[["Levine_13dim"]][which_sub_Levine_13dim, "subsampling"] <- symbol_sub
runtime_tidy[["Samusik_01"]][, "subsampling"] <- ""
runtime_tidy[["Samusik_01"]][which_sub_Samusik_01, "subsampling"] <- symbol_sub
runtime_tidy[["Samusik_all"]][, "subsampling"] <- ""
runtime_tidy[["Samusik_all"]][which_sub_Samusik_all, "subsampling"] <- symbol_sub
runtime_tidy[["Nilsson_rare"]][, "subsampling"] <- ""
runtime_tidy[["Nilsson_rare"]][which_sub_Nilsson_rare, "subsampling"] <- symbol_sub
runtime_tidy[["Mosmann_rare"]][, "subsampling"] <- ""
runtime_tidy[["Mosmann_rare"]][which_sub_Mosmann_rare, "subsampling"] <- symbol_sub

# which methods required multiple cores (from script "evaluate_runtime.R")
runtime_tidy[["Levine_32dim"]][, "cores"] <- ""
runtime_tidy[["Levine_32dim"]][which_cores_Levine_32dim, "cores"] <- symbol_cores
runtime_tidy[["Levine_13dim"]][, "cores"] <- ""
runtime_tidy[["Levine_13dim"]][which_cores_Levine_13dim, "cores"] <- symbol_cores
runtime_tidy[["Samusik_01"]][, "cores"] <- ""
runtime_tidy[["Samusik_01"]][which_cores_Samusik_01, "cores"] <- symbol_cores
runtime_tidy[["Samusik_all"]][, "cores"] <- ""
runtime_tidy[["Samusik_all"]][which_cores_Samusik_all, "cores"] <- symbol_cores
runtime_tidy[["Nilsson_rare"]][, "cores"] <- ""
runtime_tidy[["Nilsson_rare"]][which_cores_Nilsson_rare, "cores"] <- symbol_cores
runtime_tidy[["Mosmann_rare"]][, "cores"] <- ""
runtime_tidy[["Mosmann_rare"]][which_cores_Mosmann_rare, "cores"] <- symbol_cores


# bar plots: data sets with multiple populations

runtime_barplots <- vector("list", length(runtime_tidy[data_sets_multiple]))
names(runtime_barplots) <- names(runtime_tidy)[data_sets_multiple]

offsets_runtime <- list(700, 700, 850, 1000)
ymaxs_runtime <- list(37000, 35500, 43500, 52000)

x_legend <- list(4, 4.75, 4.75, 4.75)
y_legend <- list(35000, 33500, 41500, 49000)

for (i in 1:4) {
  nm <- names(runtime_tidy)[i]
  title <- paste0("Runtime: ", nm)
  filename <- paste0("../../plots/", nm, "/runtime_barplot_", nm, ".pdf")
  
  pl <- 
    ggplot(runtime_tidy[[i]], aes(x = method, y = value)) + 
    geom_bar(stat = "identity", fill = "mediumpurple") + 
    geom_text(aes(label = paste0(round(value, 0), " ", subsampling, cores), 
                  y = value + offsets_runtime[[i]], angle = 90), hjust = "left", size = 3.5) + 
    annotate("text", x = x_legend[[i]], y = y_legend[[i]], 
             label = "* subsampling required    \n^ multiple processor cores") + 
    ggtitle(title) + 
    ylim(0, ymaxs_runtime[[i]]) + 
    ylab("seconds") + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.key = element_blank(), 
          legend.title = element_blank(), 
          legend.background = element_blank())
  
  runtime_barplots[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}


# bar plots: data sets with single rare population

runtime_tidy_rare <- runtime_tidy[data_sets_single]

runtime_barplots_rare <- vector("list", length(runtime_tidy_rare))
names(runtime_barplots_rare) <- names(runtime_tidy_rare)

offsets_runtime_rare <- list(900, 1000)
ymaxs_runtime_rare <- list(43000, 49000)

x_legend_rare <- list(4.75, 4.75)
y_legend_rare <- list(40500, 46000)

for (i in 1:2) {
  nm <- names(runtime_tidy_rare)[i]
  title <- paste0("Runtime: ", nm)
  filename <- paste0("../../plots/", nm, "/runtime_barplot_", nm, ".pdf")
  
  pl <- 
    ggplot(runtime_tidy_rare[[i]], aes(x = method, y = value)) + 
    geom_bar(stat = "identity", fill = "mediumpurple") + 
    geom_text(aes(label = paste0(round(value, 0), " ", subsampling, cores), 
                  y = value + offsets_runtime_rare[[i]], angle = 90), hjust = "left", size = 3.5) + 
    annotate("text", x = x_legend_rare[[i]], y = y_legend_rare[[i]], 
             label = "* subsampling required    \n^ multiple processor cores") + 
    ggtitle(title) + 
    ylim(0, ymaxs_runtime_rare[[i]]) + 
    ylab("seconds") + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.key = element_blank(), 
          legend.title = element_blank(), 
          legend.background = element_blank())
  
  runtime_barplots_rare[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}




################################################################
### PLOTS: RUNTIME SCATTER PLOTS (RUNTIME VS. MEAN F1 SCORE) ###
################################################################

# data sets with multiple populations of interest (Levine_32dim, Levine_13dim, Samusik_01, Samusik_all)

# tidy data format (for ggplot)
f_runtime_vs_F1_tidy <- function(r, m) {
  colnames(r)[1] <- "runtime"
  r["mean_F1"] <- m[rownames(r)]
  r
}
runtime_vs_F1_tidy <- runtime_tidy[data_sets_multiple]
runtime_vs_F1_tidy <- mapply(f_runtime_vs_F1_tidy, runtime_vs_F1_tidy, mean_F1, SIMPLIFY = FALSE)

# scatter plots

runtime_scatterplots <- vector("list", length(runtime_vs_F1_tidy))
names(runtime_scatterplots) <- names(runtime_vs_F1_tidy)

ymaxs_scatter <- list(32000, 31000, 37500, 44500)

for (i in 1:4) {
  nm <- names(runtime_vs_F1_tidy)[i]
  title <- paste0("Runtime vs. mean F1: ", nm)
  filename <- paste0("../../plots/", nm, "/runtime_scatterplot_", nm, ".pdf")
  
  pl <- 
    ggplot(runtime_vs_F1_tidy[[i]], aes(x = mean_F1, y = runtime)) + 
    geom_point(shape = 4, size = 2, stroke = 1, color = "darkorchid4") + 
    geom_text_repel(aes(label = method), size = 2.5, box.padding = unit(0.3, "lines")) + 
    xlim(0, 0.8) + 
    ylim(-1500, ymaxs_scatter[[i]]) + 
    ggtitle(title) + 
    xlab("mean F1 score") + 
    ylab("runtime (seconds)") + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12))
  
  runtime_scatterplots[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, plot = pl, width = 5, height = 5)
}


# data sets with a single rare population of interest (Nilsson_rare, Mosmann_rare)

# tidy data format (for ggplot)
f_runtime_vs_F1_tidy_rare <- function(r, m) {
  colnames(r)[1] <- "runtime"
  r["F1"] <- m[rownames(r)]
  r
}
runtime_vs_F1_tidy_rare <- runtime_tidy[data_sets_single]
runtime_vs_F1_tidy_rare <- mapply(f_runtime_vs_F1_tidy_rare, runtime_vs_F1_tidy_rare, F1_rare, SIMPLIFY = FALSE)

# scatter plots

runtime_scatterplots_rare <- vector("list", length(runtime_vs_F1_tidy_rare))
names(runtime_scatterplots_rare) <- names(runtime_vs_F1_tidy_rare)

ymaxs_scatter_rare <- list(37000, 41000)
ymins_scatter_rare <- list(-3000, -3000)

for (i in 1:2) {
  nm <- names(runtime_vs_F1_tidy_rare)[i]
  title <- paste0("Runtime vs. F1: ", nm)
  filename <- paste0("../../plots/", nm, "/runtime_scatterplot_", nm, ".pdf")
  
  pl <- 
    ggplot(runtime_vs_F1_tidy_rare[[i]], aes(x = F1, y = runtime)) + 
    geom_point(shape = 4, size = 2, stroke = 1, color = "darkorchid4") + 
    geom_text_repel(aes(label = method), size = 2.5, box.padding = unit(0.3, "lines")) + 
    xlim(0, 0.8) + 
    ylim(ymins_scatter_rare[[i]], ymaxs_scatter_rare[[i]]) + 
    ggtitle(title) + 
    xlab("F1 score") + 
    ylab("runtime (seconds)") + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12))
  
  runtime_scatterplots_rare[[i]] <- pl
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
  filename <- paste0("../../plots/", nm, "/plots_multi_panel_", nm, ".pdf")
  
  pl <- ggdraw() + 
    draw_plot(barplots_mean_F1[[i]], 0.05, 0.66, 0.4, 0.33) + 
    draw_plot(boxplots_F1[[i]], 0.55, 0.66, 0.4, 0.33) + 
    draw_plot(barplots_mean_F1_pr_re[[i]], 0.05, 0.33, 0.4, 0.33) + 
    draw_plot(no_of_cells[[i]], 0.55, 0.36, 0.4, 0.30) + 
    draw_plot(runtime_barplots[[i]], 0.05, 0, 0.4, 0.33) + 
    draw_plot(runtime_scatterplots[[i]], 0.55, 0.03, 0.4, 0.30) + 
    draw_plot_label(LETTERS[1:6], 
                    c(0, 0.5, 0, 0.5, 0, 0.5), c(0.99, 0.99, 0.66, 0.66, 0.33, 0.33), size = 16)
  
  multi_panel_multiple[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, width = 13, height = 14.5)
}


# data sets with a single rare population of interest (Nilsson_rare, Mosmann_rare)

multi_panel_single <- vector("list", length(plot_data_rare))
names(multi_panel_single) <- names(plot_data_rare)

for (i in 1:2) {
  nm <- names(multi_panel_single)[i]
  filename <- paste0("../../plots/", nm, "/plots_multi_panel_", nm, ".pdf")
  
  pl <- ggdraw() + 
    draw_plot(barplots_F1_pr_re[[i]], 0.025, 0.64, 0.8, 0.35) + 
    draw_plot(runtime_barplots_rare[[i]], 0.025, 0.32, 0.8, 0.32) + 
    draw_plot(runtime_scatterplots_rare[[i]], 0.025, 0.03, 0.8, 0.29) + 
    draw_plot_label(LETTERS[1:3], c(0, 0, 0), c(0.99, 0.64, 0.32), size = 16)
  
  multi_panel_single[[i]] <- pl
  print(pl)
  
  ggplot2::ggsave(filename, pl, width = 6.5, height = 14.5)
}




###############################################
### SAVE SESSION INFORMATION AND RDATA FILE ###
###############################################

sink(file = "session_info_main_results.txt")
print(sessionInfo())
sink()

save.image("main_results.RData")



