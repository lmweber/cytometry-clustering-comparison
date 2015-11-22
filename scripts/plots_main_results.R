#########################################################################################
# R script to generate main plots and tables of results
#
# Lukas M. Weber, November 2015
#########################################################################################


library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)


# load results from previous steps

source("load_results_ACCENSE.R")
source("load_results_DensVM.R")
source("load_results_FLOCK.R")
source("load_results_PhenoGraph.R")
source("load_results_SWIFT.R")
source("load_results_truth.R")
source("load_results_all_other_methods.R")




###########################
### PREPARE DATA FRAMES ###
###########################

# combine results into lists

res_Levine <- list(ACCENSE = res_ACCENSE_Levine, 
                   DensVM = res_DensVM_Levine, 
                   FLOCK = res_FLOCK_Levine, 
                   flowMeans = res_flowMeans_Levine, 
                   FlowSOM = res_FlowSOM_Levine, 
                   FlowSOM_meta = res_FlowSOM_meta_Levine, 
                   immunoClust = res_immunoClust_Levine, 
                   immunoClust_all = res_immunoClust_all_Levine, 
                   kmeans = res_kmeans_Levine, 
                   PhenoGraph = res_PhenoGraph_Levine, 
                   Rclusterpp = res_Rclusterpp_Levine, 
                   SamSPECTRAL = res_SamSPECTRAL_Levine, 
                   SWIFT = res_SWIFT_Levine)

res_Mosmann <- list(ACCENSE = res_ACCENSE_Mosmann, 
                   DensVM = res_DensVM_Mosmann, 
                   FLOCK = res_FLOCK_Mosmann, 
                   flowMeans = res_flowMeans_Mosmann, 
                   FlowSOM = res_FlowSOM_Mosmann, 
                   FlowSOM_meta = res_FlowSOM_meta_Mosmann, 
                   immunoClust = res_immunoClust_Mosmann, 
                   immunoClust_all = res_immunoClust_all_Mosmann, 
                   kmeans = res_kmeans_Mosmann, 
                   PhenoGraph = res_PhenoGraph_Mosmann, 
                   #Rclusterpp = res_Rclusterpp_Mosmann,  # Rclusterpp did not complete for this data set
                   SamSPECTRAL = res_SamSPECTRAL_Mosmann, 
                   SWIFT = res_SWIFT_Mosmann)


# collapse into data frames

precision_df_Levine <- sapply(res_Levine, function(res) res$pr)
precision_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$pr))

recall_df_Levine <- sapply(res_Levine, function(res) res$re)
recall_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$re))

F1_df_Levine <- sapply(res_Levine, function(res) res$F1)
F1_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$F1))

labels_matched_df_Levine <- sapply(res_Levine, function(res) res$labels_matched)
labels_matched_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$labels_matched))

n_cells_df_Levine <- sapply(res_Levine, function(res) res$n_cells)
n_cells_df_Mosmann <- as.data.frame(lapply(res_Mosmann, function(res) res$n_cells))




#############################################
### MEAN F1 SCORE (LEVINE_2015_MARROW_32) ###
#############################################

n_methods <- ncol(n_cells_df_Levine)

# number of cells per true clusters

n_cells_truth_Levine <- matrix(rep(tbl_truth_Levine, n_methods), ncol = n_methods)
n_cells_truth_Levine

# mean F1 score across all true clusters, weighted by number of cells in true cluster

mean_F1_Levine <- colSums(F1_df_Levine * n_cells_truth_Levine) / colSums(n_cells_truth_Levine)
mean_F1_Levine

# arrange in descending order

ord_Levine <- rev(order(mean_F1_Levine))
mean_F1_Levine_ord <- mean_F1_Levine[ord_Levine]
mean_F1_Levine_ord

# bar plot

# tidy data format
mean_F1_Levine_tidy <- data.frame(value = mean_F1_Levine_ord)
mean_F1_Levine_tidy["method"] <- factor(rownames(mean_F1_Levine_tidy), 
                                        levels = rownames(mean_F1_Levine_tidy))
mean_F1_Levine_tidy

ggplot(mean_F1_Levine_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "blue") + 
  geom_text(aes(label = sprintf("%.2f", round(value, 2)), y = value + 0.03), size = 3.5) + 
  ylim(0, 1) + 
  xlab("") + 
  ylab("mean F1 score") + 
  ggtitle("Mean F1 score: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/results_Levine2015marrow32_mean_F1_score.pdf", height = 6, width = 5.5)




################################################
### POPULATION SIZES (LEVINE_2015_MARROW_32) ###
################################################

# plot of number of cells in each true (manually gated) population

# tidy data format
n_cells_truth_tidy <- as.data.frame(tbl_truth_Levine)
colnames(n_cells_truth_tidy) <- c("population", "value")

ggplot(n_cells_truth_tidy, aes(x = population, y = value)) + 
  geom_bar(stat = "identity", fill = "darkgray") + 
  geom_text(aes(label = value, y = value + 1300, angle = 90), size = 3.5) + 
  ylim(0, 28000) + 
  xlab("manually gated population") + 
  ylab("number of cells") + 
  ggtitle("Manually gated populations: Levine_2015_marrow_32") + 
  theme_bw()

ggsave("../plots/results_Levine2015marrow32_number_of_cells.pdf", height = 6, width = 6)




##############################################################
### HEATMAPS COMPARING ALL METHODS (LEVINE_2015_MARROW_32) ###
##############################################################

# heatmaps of F1 score, precision, and recall for each true cluster
# (with methods arranged in the same order as above)

df_Levine <- list(F1_df_Levine[, ord_Levine], 
                  precision_df_Levine[, ord_Levine], 
                  recall_df_Levine[, ord_Levine])

titles_Levine <- c("F1 score: Levine_2015_marrow_32", 
                   "Precision: Levine_2015_marrow_32", 
                   "Recall: Levine_2015_marrow_32")

filenames_Levine <- c("../plots/results_Levine2015marrow32_heatmap_F1.pdf", 
                      "../plots/results_Levine2015marrow32_heatmap_precision.pdf", 
                      "../plots/results_Levine2015marrow32_heatmap_recall.pdf")

for (i in 1:length(df_Levine)) {
  pheatmap(df_Levine[[i]], color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
           breaks = seq(0, 1, length.out = 100), display_numbers = TRUE, 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           main = titles_Levine[i], filename = filenames_Levine[i], 
           width = 6, height = 7)
}




#################################################################
### PLOTS FOR DETECTION OF RARE CELL TYPE (MOSMANN_2014_RARE) ###
#################################################################

# arrange by F1 score

ord_Mosmann <- rev(order(F1_df_Mosmann))

precision_df_Mosmann_ord <- unlist(precision_df_Mosmann[ord_Mosmann])
recall_df_Mosmann_ord <- unlist(recall_df_Mosmann[ord_Mosmann])
F1_df_Mosmann_ord <- unlist(F1_df_Mosmann[ord_Mosmann])

# tidy data format

plot_data_Mosmann <- data.frame(precision = precision_df_Mosmann_ord, 
                                recall = recall_df_Mosmann_ord, 
                                F1_score = F1_df_Mosmann_ord)
plot_data_Mosmann["method"] <- factor(rownames(plot_data_Mosmann), 
                                      levels = rownames(plot_data_Mosmann))

plot_data_Mosmann <- melt(plot_data_Mosmann, 
                          id.vars = "method", 
                          measure.vars = c("F1_score", "precision", "recall"))

# bar plots

ggplot(plot_data_Mosmann, aes(x = method, y = value, group = variable, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0, 1) + 
  xlab("") + 
  ylab("") + 
  ggtitle("Rare cell population: Mosmann_2014_rare") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.title = element_blank())

ggsave("../plots/results_Mosmann2014rare_F1_precision_recall.pdf", height = 6, width = 7)


