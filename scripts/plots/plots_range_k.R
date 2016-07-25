#########################################################################################
# R script to generate plot of results for FlowSOM_meta using range of values for number
# of clusters k
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)
library(ggplot2)
library(reshape2)

# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")

# color-blind friendly palettes (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cb_pal_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_pal_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




####################
### LOAD RESULTS ###
####################

RES_DIR <- "../results_range_k"

# load results

file_res_Levine_32 <- file.path(RES_DIR, "results_FlowSOM_meta_range_k_Levine_32.txt")
file_res_Levine_13 <- file.path(RES_DIR, "results_FlowSOM_meta_range_k_Levine_13.txt")
file_res_Nilsson <- file.path(RES_DIR, "results_FlowSOM_meta_range_k_Nilsson.txt")
file_res_Mosmann <- file.path(RES_DIR, "results_FlowSOM_meta_range_k_Mosmann.txt")

res_Levine_32 <- read.table(file_res_Levine_32, header = TRUE, sep = "\t", comment.char = "")
res_Levine_13 <- read.table(file_res_Levine_13, header = TRUE, sep = "\t", comment.char = "")
res_Nilsson <- read.table(file_res_Nilsson, header = TRUE, sep = "\t", comment.char = "")
res_Mosmann <- read.table(file_res_Mosmann, header = TRUE, sep = "\t", comment.char = "")

res_Levine_32
res_Levine_13
res_Nilsson
res_Mosmann




############
### PLOT ###
############

# plot mean F1 score or F1 score vs. k

# tidy data format for ggplot

res_F1 <- data.frame(k = res_Levine_32$k, 
                     Levine_2015_marrow_32 = res_Levine_32$mean_F1, 
                     Levine_2015_marrow_13 = res_Levine_13$mean_F1, 
                     Nilsson_2013_HSC = res_Nilsson$F1, 
                     Mosmann_2014_activ = res_Mosmann$F1)

res_F1_tidy <- melt(res_F1, id.vars = "k", measure_vars = c("Levine_32", "Levine_13", "Nilsson", "Mosmann"))


# line plot

lineplot_F1_range_k <- 
  ggplot(res_F1_tidy, aes(x = k, y = value, group = variable)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(color = variable)) + 
  geom_vline(xintercept = 40, color = "red", linetype = 1) + 
  scale_color_manual(values = cb_pal_black[1:4]) + 
  ylim(0, 0.82) + 
  scale_x_continuous(breaks = seq(5, 80, by = 5)) + 
  ylab("mean F1 score or F1 score") + 
  xlab("number of clusters k") + 
  ggtitle("FlowSOM_meta") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        legend.position = c(0.80, 0.15), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_rect("white"))

lineplot_F1_range_k
ggplot2::ggsave("../plots/plots_range_k/results_FlowSOM_meta_F1_range_k.pdf", width = 6, height = 5)



