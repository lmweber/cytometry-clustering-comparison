#########################################################################################
# R script to generate plots for stability analysis
#
# Lukas M. Weber, February 2016
#########################################################################################


library(ggplot2)
library(reshape2)




#################
### LOAD DATA ###
#################

STABILITY_DIR <- "../results/stability_analysis"


# read saved data

files_FLOCK <- list.files(file.path(STABILITY_DIR, "FLOCK"), "^stability", full.names = TRUE)
files_FLOCK <- files_FLOCK[c(2, 1, 4, 3)]
res_stability_FLOCK <- lapply(files_FLOCK, function(f) read.table(f, header = TRUE, sep = "\t"))

files_flowMeans <- list.files(file.path(STABILITY_DIR, "flowMeans"), "^stability", full.names = TRUE)
files_flowMeans <- files_flowMeans[c(2, 1, 4, 3)]
res_stability_flowMeans <- lapply(files_flowMeans, function(f) read.table(f, header = TRUE, sep = "\t"))

files_FlowSOM <- list.files(file.path(STABILITY_DIR, "FlowSOM"), "^stability", full.names = TRUE)
files_FlowSOM <- files_FlowSOM[c(2, 1, 4, 3)]
res_stability_FlowSOM <- lapply(files_FlowSOM, function(f) read.table(f, header = TRUE, sep = "\t"))

files_FlowSOM_meta <- list.files(file.path(STABILITY_DIR, "FlowSOM_meta"), "^stability", full.names = TRUE)
files_FlowSOM_meta <- files_FlowSOM_meta[c(2, 1, 4, 3)]
res_stability_FlowSOM_meta <- lapply(files_FlowSOM_meta, function(f) read.table(f, header = TRUE, sep = "\t"))

files_immunoClust <- list.files(file.path(STABILITY_DIR, "immunoClust"), "^stability", full.names = TRUE)
files_immunoClust <- files_immunoClust[c(2, 1, 4, 3)]
res_stability_immunoClust <- lapply(files_immunoClust, function(f) read.table(f, header = TRUE, sep = "\t"))

files_immunoClust_all <- list.files(file.path(STABILITY_DIR, "immunoClust_all"), "^stability", full.names = TRUE)
files_immunoClust_all <- files_immunoClust_all[c(2, 1, 4, 3)]
res_stability_immunoClust_all <- lapply(files_immunoClust_all, function(f) read.table(f, header = TRUE, sep = "\t"))

files_PhenoGraph <- list.files(file.path(STABILITY_DIR, "PhenoGraph"), "^stability", full.names = TRUE)
files_PhenoGraph <- files_PhenoGraph[c(2, 1, 4, 3)]
res_stability_PhenoGraph <- lapply(files_PhenoGraph, function(f) read.table(f, header = TRUE, sep = "\t"))

files_SamSPECTRAL <- list.files(file.path(STABILITY_DIR, "SamSPECTRAL_without_Mosmann"), "^stability", full.names = TRUE)
files_SamSPECTRAL <- files_SamSPECTRAL[c(2, 1, 3)]  # Mosmann not complete yet
res_stability_SamSPECTRAL <- lapply(files_SamSPECTRAL, function(f) read.table(f, header = TRUE, sep = "\t"))


# combine results into lists

res_stability_Levine_32 <- list(FLOCK = res_stability_FLOCK[[1]], 
                                flowMeans = res_stability_flowMeans[[1]], 
                                FlowSOM = res_stability_FlowSOM[[1]], 
                                FlowSOM_meta = res_stability_FlowSOM_meta[[1]], 
                                immunoClust = res_stability_immunoClust[[1]], 
                                immunoClust_all = res_stability_immunoClust_all[[1]], 
                                PhenoGraph = res_stability_PhenoGraph[[1]], 
                                SamSPECTRAL = res_stability_SamSPECTRAL[[1]])

res_stability_Levine_13 <- list(FLOCK = res_stability_FLOCK[[2]], 
                                flowMeans = res_stability_flowMeans[[2]], 
                                FlowSOM = res_stability_FlowSOM[[2]], 
                                FlowSOM_meta = res_stability_FlowSOM_meta[[2]], 
                                immunoClust = res_stability_immunoClust[[2]], 
                                immunoClust_all = res_stability_immunoClust_all[[2]], 
                                PhenoGraph = res_stability_PhenoGraph[[2]], 
                                SamSPECTRAL = res_stability_SamSPECTRAL[[2]])

res_stability_Nilsson <- list(FLOCK = res_stability_FLOCK[[3]], 
                              flowMeans = res_stability_flowMeans[[3]], 
                              FlowSOM = res_stability_FlowSOM[[3]], 
                              FlowSOM_meta = res_stability_FlowSOM_meta[[3]], 
                              immunoClust = res_stability_immunoClust[[3]], 
                              immunoClust_all = res_stability_immunoClust_all[[3]], 
                              PhenoGraph = res_stability_PhenoGraph[[3]], 
                              SamSPECTRAL = res_stability_SamSPECTRAL[[3]])

res_stability_Mosmann <- list(FLOCK = res_stability_FLOCK[[4]], 
                              flowMeans = res_stability_flowMeans[[4]], 
                              FlowSOM = res_stability_FlowSOM[[4]], 
                              FlowSOM_meta = res_stability_FlowSOM_meta[[4]], 
                              immunoClust = res_stability_immunoClust[[4]], 
                              immunoClust_all = res_stability_immunoClust_all[[4]], 
                              PhenoGraph = res_stability_PhenoGraph[[4]]) #, 
                              #SamSPECTRAL = res_stability_SamSPECTRAL[[4]])  # not complete yet


# collapse into data frames

df_stability_F1_Levine_32 <- as.data.frame(sapply(res_stability_Levine_32, function(res) res$mean_F1))
df_stability_F1_Levine_13 <- as.data.frame(sapply(res_stability_Levine_13, function(res) res$mean_F1))
df_stability_F1_Nilsson <- as.data.frame(sapply(res_stability_Nilsson, function(res) res$F1))
df_stability_F1_Mosmann <- as.data.frame(sapply(res_stability_Mosmann, function(res) res$F1))

df_stability_pr_Levine_32 <- as.data.frame(sapply(res_stability_Levine_32, function(res) res$mean_pr))
df_stability_pr_Levine_13 <- as.data.frame(sapply(res_stability_Levine_13, function(res) res$mean_pr))
df_stability_pr_Nilsson <- as.data.frame(sapply(res_stability_Nilsson, function(res) res$pr))
df_stability_pr_Mosmann <- as.data.frame(sapply(res_stability_Mosmann, function(res) res$pr))

df_stability_re_Levine_32 <- as.data.frame(sapply(res_stability_Levine_32, function(res) res$mean_re))
df_stability_re_Levine_13 <- as.data.frame(sapply(res_stability_Levine_13, function(res) res$mean_re))
df_stability_re_Nilsson <- as.data.frame(sapply(res_stability_Nilsson, function(res) res$re))
df_stability_re_Mosmann <- as.data.frame(sapply(res_stability_Mosmann, function(res) res$re))

df_stability_runtime_Levine_32 <- as.data.frame(sapply(res_stability_Levine_32, function(res) res$runtime))
df_stability_runtime_Levine_13 <- as.data.frame(sapply(res_stability_Levine_13, function(res) res$runtime))
df_stability_runtime_Nilsson <- as.data.frame(sapply(res_stability_Nilsson, function(res) res$runtime))
df_stability_runtime_Mosmann <- as.data.frame(sapply(res_stability_Mosmann, function(res) res$runtime))


# arrange columns in same order as previously (excluding methods with no stability analysis)

ord_stability_Levine_32 <- c(4, 2, 1, 7, 3, 8, 6, 5)
ord_stability_Levine_13 <- c(4, 7, 2, 1, 8, 3, 6, 5)
ord_stability_Nilsson <- c(2, 4, 3, 7, 6, 1, 8, 5)
ord_stability_Mosmann <- c(7, 4, 3, 2, 1, 6, 5)  # SamSPECTRAL not included yet

df_stability_F1_Levine_32 <- df_stability_F1_Levine_32[, ord_stability_Levine_32]
df_stability_F1_Levine_13 <- df_stability_F1_Levine_13[, ord_stability_Levine_13]
df_stability_F1_Nilsson <- df_stability_F1_Nilsson[, ord_stability_Nilsson]
df_stability_F1_Mosmann <- df_stability_F1_Mosmann[, ord_stability_Mosmann]

df_stability_pr_Levine_32 <- df_stability_pr_Levine_32[, ord_stability_Levine_32]
df_stability_pr_Levine_13 <- df_stability_pr_Levine_13[, ord_stability_Levine_13]
df_stability_pr_Nilsson <- df_stability_pr_Nilsson[, ord_stability_Nilsson]
df_stability_pr_Mosmann <- df_stability_pr_Mosmann[, ord_stability_Mosmann]

df_stability_re_Levine_32 <- df_stability_re_Levine_32[, ord_stability_Levine_32]
df_stability_re_Levine_13 <- df_stability_re_Levine_13[, ord_stability_Levine_13]
df_stability_re_Nilsson <- df_stability_re_Nilsson[, ord_stability_Nilsson]
df_stability_re_Mosmann <- df_stability_re_Mosmann[, ord_stability_Mosmann]

df_stability_runtime_Levine_32 <- df_stability_runtime_Levine_32[, ord_stability_Levine_32]
df_stability_runtime_Levine_13 <- df_stability_runtime_Levine_13[, ord_stability_Levine_13]
df_stability_runtime_Nilsson <- df_stability_runtime_Nilsson[, ord_stability_Nilsson]
df_stability_runtime_Mosmann <- df_stability_runtime_Mosmann[, ord_stability_Mosmann]




##############################################
### TIDY DATA: F1 SCORE, PRECISION, RECALL ###
##############################################

# arrange in tidy data format for ggplot2


# Levine_32
df_stability_Levine_32_tidy <- data.frame(F1_score = as.vector(as.matrix(df_stability_F1_Levine_32)), 
                                          precision = as.vector(as.matrix(df_stability_pr_Levine_32)), 
                                          recall = as.vector(as.matrix(df_stability_re_Levine_32)))
df_stability_Levine_32_tidy["method"] <- rep(factor(colnames(df_stability_F1_Levine_32), 
                                                    levels = colnames(df_stability_F1_Levine_32)), 
                                             each = nrow(df_stability_F1_Levine_32))
df_stability_Levine_32_tidy <- melt(df_stability_Levine_32_tidy, 
                                    id.vars = "method", 
                                    measure.vars = c("F1_score", "precision", "recall"))

# Levine_13
df_stability_Levine_13_tidy <- data.frame(F1_score = as.vector(as.matrix(df_stability_F1_Levine_13)), 
                                          precision = as.vector(as.matrix(df_stability_pr_Levine_13)), 
                                          recall = as.vector(as.matrix(df_stability_re_Levine_13)))
df_stability_Levine_13_tidy["method"] <- rep(factor(colnames(df_stability_F1_Levine_13), 
                                                    levels = colnames(df_stability_F1_Levine_13)), 
                                             each = nrow(df_stability_F1_Levine_13))
df_stability_Levine_13_tidy <- melt(df_stability_Levine_13_tidy, 
                                    id.vars = "method", 
                                    measure.vars = c("F1_score", "precision", "recall"))

# Nilsson
df_stability_Nilsson_tidy <- data.frame(F1_score = as.vector(as.matrix(df_stability_F1_Nilsson)), 
                                        precision = as.vector(as.matrix(df_stability_pr_Nilsson)), 
                                        recall = as.vector(as.matrix(df_stability_re_Nilsson)))
df_stability_Nilsson_tidy["method"] <- rep(factor(colnames(df_stability_F1_Nilsson), 
                                                  levels = colnames(df_stability_F1_Nilsson)), 
                                           each = nrow(df_stability_F1_Nilsson))
df_stability_Nilsson_tidy <- melt(df_stability_Nilsson_tidy, 
                                  id.vars = "method", 
                                  measure.vars = c("F1_score", "precision", "recall"))

# Mosmann
df_stability_Mosmann_tidy <- data.frame(F1_score = as.vector(as.matrix(df_stability_F1_Mosmann)), 
                                        precision = as.vector(as.matrix(df_stability_pr_Mosmann)), 
                                        recall = as.vector(as.matrix(df_stability_re_Mosmann)))
df_stability_Mosmann_tidy["method"] <- rep(factor(colnames(df_stability_F1_Mosmann), 
                                                  levels = colnames(df_stability_F1_Mosmann)), 
                                           each = nrow(df_stability_F1_Mosmann))
df_stability_Mosmann_tidy <- melt(df_stability_Mosmann_tidy, 
                                  id.vars = "method", 
                                  measure.vars = c("F1_score", "precision", "recall"))




##########################
### TIDY DATA: RUNTIME ###
##########################

# arrange in tidy data format for ggplot2

# Levine_32: runtime
df_stability_runtime_Levine_32_tidy <- data.frame(runtime = as.vector(as.matrix(df_stability_runtime_Levine_32)))
df_stability_runtime_Levine_32_tidy["method"] <- rep(factor(colnames(df_stability_runtime_Levine_32), 
                                                            levels = colnames(df_stability_runtime_Levine_32)), 
                                                     each = nrow(df_stability_runtime_Levine_32))
df_stability_runtime_Levine_32_tidy <- melt(df_stability_runtime_Levine_32_tidy, 
                                            id.vars = "method", measure.vars = "runtime")

# Levine_13: runtime
df_stability_runtime_Levine_13_tidy <- data.frame(runtime = as.vector(as.matrix(df_stability_runtime_Levine_13)))
df_stability_runtime_Levine_13_tidy["method"] <- rep(factor(colnames(df_stability_runtime_Levine_13), 
                                                            levels = colnames(df_stability_runtime_Levine_13)), 
                                                     each = nrow(df_stability_runtime_Levine_13))
df_stability_runtime_Levine_13_tidy <- melt(df_stability_runtime_Levine_13_tidy, 
                                            id.vars = "method", measure.vars = "runtime")

# Nilsson: runtime
df_stability_runtime_Nilsson_tidy <- data.frame(runtime = as.vector(as.matrix(df_stability_runtime_Nilsson)))
df_stability_runtime_Nilsson_tidy["method"] <- rep(factor(colnames(df_stability_runtime_Nilsson), 
                                                          levels = colnames(df_stability_runtime_Nilsson)), 
                                                   each = nrow(df_stability_runtime_Nilsson))
df_stability_runtime_Nilsson_tidy <- melt(df_stability_runtime_Nilsson_tidy, 
                                          id.vars = "method", measure.vars = "runtime")

# Mosmann: runtime
df_stability_runtime_Mosmann_tidy <- data.frame(runtime = as.vector(as.matrix(df_stability_runtime_Mosmann)))
df_stability_runtime_Mosmann_tidy["method"] <- rep(factor(colnames(df_stability_runtime_Mosmann), 
                                                          levels = colnames(df_stability_runtime_Mosmann)), 
                                                   each = nrow(df_stability_runtime_Mosmann))
df_stability_runtime_Mosmann_tidy <- melt(df_stability_runtime_Mosmann_tidy, 
                                          id.vars = "method", measure.vars = "runtime")


# remove PhenoGraph (not comparable since used Python implementation and multiple cores)
# note: ggplot2 doesn't require fixing unused factor levels or re-numbering rows

df_stability_runtime_Levine_32_tidy <- 
  df_stability_runtime_Levine_32_tidy[df_stability_runtime_Levine_32_tidy$method != "PhenoGraph", ]

df_stability_runtime_Levine_13_tidy <- 
  df_stability_runtime_Levine_13_tidy[df_stability_runtime_Levine_13_tidy$method != "PhenoGraph", ]

df_stability_runtime_Nilsson_tidy <- 
  df_stability_runtime_Nilsson_tidy[df_stability_runtime_Nilsson_tidy$method != "PhenoGraph", ]

df_stability_runtime_Mosmann_tidy <- 
  df_stability_runtime_Mosmann_tidy[df_stability_runtime_Mosmann_tidy$method != "PhenoGraph", ]




##############################################
### BOX PLOTS: F1 SCORE, PRECISION, RECALL ###
##############################################

# box plots of mean F1 score, mean precision, mean recall
# (for data sets with multiple populations of interest: Levine_2015_marrow_32, Levine_2015_marrow_13)


# Levine_32
boxplots_stability_Levine_32 <- 
  ggplot(df_stability_Levine_32_tidy, aes(x = method, y = value, color = variable, fill = variable)) + 
  geom_boxplot(width = 0.85, position = position_dodge(0.75), outlier.size = 0.75) + 
  scale_fill_hue(c = 60, l = 90) + 
  ylim(0, 1) + 
  ggtitle("Stability of clustering results: Levine_2015_marrow_32") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.13, 0.15), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

boxplots_stability_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/stability_analysis/stability_boxplots_Levine2015marrow32.pdf", 
                width = 5, height = 5)


# Levine_13
boxplots_stability_Levine_13 <- 
  ggplot(df_stability_Levine_13_tidy, aes(x = method, y = value, color = variable, fill = variable)) + 
  geom_boxplot(width = 0.85, position = position_dodge(0.75), outlier.size = 0.75) + 
  scale_fill_hue(c = 60, l = 90) + 
  ylim(0, 1) + 
  ggtitle("Stability of clustering results: Levine_2015_marrow_13") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.13, 0.15), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

boxplots_stability_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/stability_analysis/stability_boxplots_Levine2015marrow13.pdf", 
                width = 5, height = 5)


# box plots of F1 score, precision, recall
# (for data sets with a single rare cell population of interest: Nilsson_2013_HSC, Mosmann_2014_activ)


# Nilsson
boxplots_stability_Nilsson <- 
  ggplot(df_stability_Nilsson_tidy, aes(x = method, y = value, color = variable, fill = variable)) + 
  geom_boxplot(width = 0.85, position = position_dodge(0.75), outlier.size = 0.75) + 
  scale_fill_hue(c = 60, l = 90) + 
  ylim(0, 1) + 
  ggtitle("Stability of clustering results: Nilsson_2013_HSC") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.12, 0.12), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

boxplots_stability_Nilsson
ggplot2::ggsave("../plots/Nilsson_2013_HSC/stability_analysis/stability_boxplots_Nilsson2013HSC.pdf", 
                width = 5, height = 5)


# Mosmann
boxplots_stability_Mosmann <- 
  ggplot(df_stability_Mosmann_tidy, aes(x = method, y = value, color = variable, fill = variable)) + 
  geom_boxplot(width = 0.85, position = position_dodge(0.75), outlier.size = 0.75) + 
  scale_fill_hue(c = 60, l = 90) + 
  ylim(0, 1) + 
  ggtitle("Stability of clustering results: Mosmann_2014_activ") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.89, 0.9), 
        legend.key.size = unit(5, "mm"), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank())

boxplots_stability_Mosmann
ggplot2::ggsave("../plots/Mosmann_2014_activ/stability_analysis/stability_boxplots_Mosmann2014activ.pdf", 
                width = 5, height = 5)




##########################
### BOX PLOTS: RUNTIME ###
##########################

# box plots of runtime

# Levine_32
boxplots_stability_runtime_Levine_32 <- 
  ggplot(df_stability_runtime_Levine_32_tidy, aes(x = method, y = value)) + 
  geom_boxplot(width = 0.4, outlier.size = 0.8, color = "purple4", fill = "mediumpurple1") + 
  ggtitle("Distribution of runtimes: Levine_2015_marrow_32") + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

boxplots_stability_runtime_Levine_32
ggplot2::ggsave("../plots/Levine_2015_marrow_32/stability_analysis/stability_boxplots_runtime_Levine2015marrow32.pdf", 
                width = 5, height = 5)


# Levine_13
boxplots_stability_runtime_Levine_13 <- 
  ggplot(df_stability_runtime_Levine_13_tidy, aes(x = method, y = value)) + 
  geom_boxplot(width = 0.4, outlier.size = 0.8, color = "purple4", fill = "mediumpurple1") + 
  ggtitle("Distribution of runtimes: Levine_2015_marrow_13") + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

boxplots_stability_runtime_Levine_13
ggplot2::ggsave("../plots/Levine_2015_marrow_13/stability_analysis/stability_boxplots_runtime_Levine2015marrow13.pdf", 
                width = 5, height = 5)


# Nilsson
boxplots_stability_runtime_Nilsson <- 
  ggplot(df_stability_runtime_Nilsson_tidy, aes(x = method, y = value)) + 
  geom_boxplot(width = 0.4, outlier.size = 0.8, color = "purple4", fill = "mediumpurple1") + 
  ggtitle("Distribution of runtimes: Nilsson_2013_HSC") + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

boxplots_stability_runtime_Nilsson
ggplot2::ggsave("../plots/Nilsson_2013_HSC/stability_analysis/stability_boxplots_runtime_Nilsson2013HSC.pdf", 
                width = 5, height = 5)


# Mosmann
boxplots_stability_runtime_Mosmann <- 
  ggplot(df_stability_runtime_Mosmann_tidy, aes(x = method, y = value)) + 
  geom_boxplot(width = 0.4, outlier.size = 0.8, color = "purple4", fill = "mediumpurple1") + 
  ggtitle("Distribution of runtimes: Mosmann_2014_activ") + 
  ylab("seconds") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

boxplots_stability_runtime_Mosmann
ggplot2::ggsave("../plots/Mosmann_2014_activ/stability_analysis/stability_boxplots_runtime_Mosmann2014activ.pdf", 
                width = 5, height = 5)



