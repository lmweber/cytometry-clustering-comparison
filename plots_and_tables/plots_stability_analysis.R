#########################################################################################
# R script to generate plots for stability analysis
#
# Lukas Weber, September 2016
#########################################################################################


library(ggplot2)
library(reshape2)
library(cowplot)  # note masks ggplot2::ggsave()
library(colorspace)

# color-blind friendly palettes (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cb_pal_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_pal_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ggplot2 default palette (access with ggplot_build(p)$data)
gg_pal <- c("#F8766D", "#00BA38", "#619CFF")

# custom palette: adjust chroma, luminance
custom_pal <- c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])
hcl <- coords(as(hex2RGB(custom_pal), "polarLUV"))
hcl[, "L"] <- 85
hcl[, "C"] <- 40
custom_pal <- hex(polarLUV(hcl), fixup = TRUE)




#################
### LOAD DATA ###
#################

DIR_RANDOM_STARTS <- "../../results/stability/random_starts"
DIR_BOOTSTRAP <- "../../results/stability/bootstrap"

datasets <- c("Levine_32dim", "Mosmann_rare")

files_random_starts <- list(
  FLOCK = paste0(DIR_RANDOM_STARTS, "/FLOCK/stability_random_starts_FLOCK_", datasets, ".txt"), 
  flowMeans = paste0(DIR_RANDOM_STARTS, "/flowMeans/stability_random_starts_flowMeans_", datasets, ".txt"), 
  flowPeaks = paste0(DIR_RANDOM_STARTS, "/flowPeaks/stability_random_starts_flowPeaks_", datasets, ".txt"), 
  FlowSOM = paste0(DIR_RANDOM_STARTS, "/FlowSOM/stability_random_starts_flowSOM_", datasets, ".txt"), 
  FlowSOM_pre = paste0(DIR_RANDOM_STARTS, "/FlowSOM_pre/stability_random_starts_flowSOM_pre_", datasets, ".txt"), 
  kmeans = paste0(DIR_RANDOM_STARTS, "/kmeans/stability_random_starts_kmeans_", datasets, ".txt"), 
  SamSPECTRAL = paste0(DIR_RANDOM_STARTS, "/SamSPECTRAL/stability_random_starts_SamSPECTRAL_", datasets, ".txt")
)

files_bootstrap <- list(
  FLOCK = paste0(DIR_BOOTSTRAP, "/FLOCK/stability_bootstrap_FLOCK_", datasets, ".txt"), 
  flowMeans = paste0(DIR_BOOTSTRAP, "/flowMeans/stability_bootstrap_flowMeans_", datasets, ".txt"), 
  flowPeaks = paste0(DIR_BOOTSTRAP, "/flowPeaks/stability_bootstrap_flowPeaks_", datasets, ".txt"), 
  FlowSOM = paste0(DIR_BOOTSTRAP, "/FlowSOM/stability_bootstrap_FlowSOM_", datasets, ".txt"), 
  FlowSOM_pre = paste0(DIR_BOOTSTRAP, "/FlowSOM_pre/stability_bootstrap_FlowSOM_pre_", datasets, ".txt"), 
  kmeans = paste0(DIR_BOOTSTRAP, "/kmeans/stability_bootstrap_kmeans_", datasets, ".txt"), 
  SamSPECTRAL = paste0(DIR_BOOTSTRAP, "/SamSPECTRAL/stability_bootstrap_SamSPECTRAL_", datasets, ".txt")
)


# rearrange: data sets at top level of list

fn_rearrange_files <- function(files, datasets) {
  methods <- names(files)
  tmp <- vector("list", length(datasets))
  names(tmp) <- datasets
  for (i in 1:length(datasets)) {
    tmp[[i]] <- vector("list", length(methods))
    names(tmp[[i]]) <- methods
    for (j in 1:length(files)) {
      tmp[[i]][[j]] <- files[[j]][i]
    }
  }
  tmp
}

files_random_starts <- fn_rearrange_files(files_random_starts, datasets)
files_bootstrap <- fn_rearrange_files(files_bootstrap, datasets)


# load data

data_random_starts <- lapply(files_random_starts, function(d) {
  lapply(d, read.table, header = TRUE, sep = "\t")
})

data_bootstrap <- lapply(files_bootstrap, function(d) {
  lapply(d, read.table, header = TRUE, sep = "\t")
})




###########################
### PREPARE DATA FRAMES ###
###########################

# collapse into data frames for each of F1, pr, re

fn_rearrange_data <- function(data) {
  datasets <- names(data)
  tmp <- vector("list", length(data))
  names(tmp) <- datasets
  for (i in 1:length(tmp)) {
    tmp[[i]] <- vector("list", 3)
    names(tmp[[i]]) <- c("pr", "re", "F1")  ## note ordering: pr, re, F1 (instead of F1 first)
    for (j in 1:length(tmp[[i]])) {
      tmp[[i]][[j]] <- sapply(data[[i]], function(d) d[, j])
    }
  }
  tmp
}

df_random_starts <- fn_rearrange_data(data_random_starts)
df_bootstrap <- fn_rearrange_data(data_bootstrap)


# separate by data set; arrange clustering methods in same order as previously

ord_stability_Levine_32dim <- c("FlowSOM", "flowMeans", "FLOCK", "SamSPECTRAL", "FlowSOM_pre", "kmeans", "flowPeaks")
ord_stability_Mosmann_rare <- c("FlowSOM_pre", "FlowSOM", "SamSPECTRAL", "flowMeans", "kmeans", "FLOCK", "flowPeaks")

df_random_starts_Levine_32dim <- lapply(df_random_starts[[1]], function(d) d[, ord_stability_Levine_32dim])[c(3, 1, 2)]
df_random_starts_Mosmann_rare <- lapply(df_random_starts[[2]], function(d) d[, ord_stability_Mosmann_rare])[c(3, 1, 2)]

df_bootstrap_Levine_32dim <- lapply(df_bootstrap[[1]], function(d) d[, ord_stability_Levine_32dim])[c(3, 1, 2)]
df_bootstrap_Mosmann_rare <- lapply(df_bootstrap[[2]], function(d) d[, ord_stability_Mosmann_rare])[c(3, 1, 2)]


# tidy data format (for ggplot)

f_plot_data <- function(df) {
  d <- data.frame(F1_score = as.vector(as.matrix(df[["F1"]])), 
                  precision = as.vector(as.matrix(df[["pr"]])), 
                  recall = as.vector(as.matrix(df[["re"]])))
  d["method"] <- rep(factor(colnames(df[[1]]), levels = colnames(df[[1]])), each = nrow(df[[1]]))
  d <- melt(d, id.vars = "method", measure.vars = c("F1_score", "precision", "recall"))
  d
}

plot_data_random_starts_Levine_32dim <- f_plot_data(df_random_starts_Levine_32dim)
plot_data_random_starts_Mosmann_rare <- f_plot_data(df_random_starts_Mosmann_rare)

plot_data_bootstrap_Levine_32dim <- f_plot_data(df_bootstrap_Levine_32dim)
plot_data_bootstrap_Mosmann_rare <- f_plot_data(df_bootstrap_Mosmann_rare)




#####################################
### BOX PLOTS: STABILITY ANALYSIS ###
#####################################

# box plots for stability analysis (both random starts and bootstrap)


plot_data <- list(
  plot_data_random_starts_Levine_32dim, 
  plot_data_random_starts_Mosmann_rare, 
  plot_data_bootstrap_Levine_32dim, 
  plot_data_bootstrap_Mosmann_rare
)

titles <- c(
  "Stability (random starts): Levine_32dim", 
  "Stability (random starts): Mosmann_rare", 
  "Stability (bootstrap): Levine_32dim", 
  "Stability (bootstrap): Mosmann_rare"
)

filenames <- c(
  "../../plots/Levine_32dim/stability_analysis/stability_random_starts_Levine_32dim.pdf", 
  "../../plots/Mosmann_rare/stability_analysis/stability_random_starts_Mosmann_rare.pdf", 
  "../../plots/Levine_32dim/stability_analysis/stability_bootstrap_Levine_32dim.pdf", 
  "../../plots/Mosmann_rare/stability_analysis/stability_bootstrap_Mosmann_rare.pdf"
)


for (i in 1:4) {
  pl <- 
    ggplot(plot_data[[i]], 
           aes(x = method, y = value, color = variable, fill = variable)) + 
    geom_boxplot(width = 0.85, position = position_dodge(0.75), outlier.size = 0.75) + 
    scale_fill_manual(values = custom_pal) + 
    scale_color_manual(values = c(gg_pal[1], cb_pal_black[4], cb_pal_black[3])) + 
    ylim(0, 1) + 
    ggtitle(titles[i]) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "right", 
          legend.key.size = unit(5, "mm"), 
          legend.key = element_blank(), 
          legend.title = element_blank(), 
          legend.background = element_blank())
  
  print(pl)
  
  ggplot2::ggsave(filenames[i], plot = pl, width = 6, height = 5)
}



