#########################################################################################
# R script to run stability analysis for methods in parallel
#
# This script runs clustering methods multiple times in parallel using BiocParallel.
# Additional methods are run in series (since they do not return results as R objects) in
# "stability_analysis_series.R".
# 
# Note: Set random seeds for each bplapply call, but not for individual iterations.
# 
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)
library(BiocParallel)

# helper functions
source("../helper_match_clusters_and_evaluate.R")
source("../helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels for stability analysis
source("load_results_truth_stability.R")

# number of times to run each method
n <- 30

# directory to save results
RESULTS_DIR <- "../../results_stability_analysis"




#################
### LOAD DATA ###
#################

DATA_DIR <- "../../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")

data_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_Levine_32, transformation = FALSE))
data_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_Levine_13, transformation = FALSE))
data_Nilsson <- flowCore::exprs(flowCore::read.FCS(file_Nilsson, transformation = FALSE))
data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)


# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- c(7:9, 11:21)

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)


# subset data

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
data_Nilsson <- data_Nilsson[, marker_cols_Nilsson]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)




#################
### flowMeans ###
#################

library(flowMeans)

source("run_flowMeans_stability.R")

# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32, data_Levine_13, data_Nilsson, data_Mosmann)), n)

res_stability_flowMeans <- bplapply(data_list, run_flowMeans_stability, 
                                    BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_flowMeans, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_flowMeans, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_flowMeans, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_flowMeans, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "flowMeans", "stability_flowMeans_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




###############
### FlowSOM ###
###############

library(FlowSOM)

source("run_FlowSOM_stability.R")


# FlowSOM requires input data as flowFrame objects

data_Levine_32_FlowSOM <- flowCore::read.FCS(file_Levine_32, transformation = FALSE)
data_Levine_13_FlowSOM <- flowCore::read.FCS(file_Levine_13, transformation = FALSE)
data_Nilsson_FlowSOM <- flowCore::read.FCS(file_Nilsson, transformation = FALSE)
data_Mosmann_FlowSOM <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)


# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32_FlowSOM, data_Levine_13_FlowSOM, 
                           data_Nilsson_FlowSOM, data_Mosmann_FlowSOM)), n)

res_stability_FlowSOM <- bplapply(data_list, run_FlowSOM_stability, 
                                  BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_FlowSOM, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_FlowSOM, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_FlowSOM, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_FlowSOM, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "FlowSOM", "stability_FlowSOM_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




####################
### FlowSOM_meta ###
####################

library(FlowSOM)

source("run_FlowSOM_meta_stability.R")


# FlowSOM_meta requires input data as flowFrame objects

data_Levine_32_FlowSOM_meta <- flowCore::read.FCS(file_Levine_32, transformation = FALSE)
data_Levine_13_FlowSOM_meta <- flowCore::read.FCS(file_Levine_13, transformation = FALSE)
data_Nilsson_FlowSOM_meta <- flowCore::read.FCS(file_Nilsson, transformation = FALSE)
data_Mosmann_FlowSOM_meta <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)


# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32_FlowSOM_meta, data_Levine_13_FlowSOM_meta, 
                           data_Nilsson_FlowSOM_meta, data_Mosmann_FlowSOM_meta)), n)

res_stability_FlowSOM_meta <- bplapply(data_list, run_FlowSOM_meta_stability, 
                                       BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_FlowSOM_meta, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_FlowSOM_meta, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_FlowSOM_meta, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_FlowSOM_meta, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "FlowSOM_meta", "stability_FlowSOM_meta_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




###############
### k-means ###
###############

source("run_kmeans_stability.R")

# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32, data_Levine_13, data_Nilsson, data_Mosmann)), n)

res_stability_kmeans <- bplapply(data_list, run_kmeans_stability, 
                                 BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_kmeans, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_kmeans, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_kmeans, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_kmeans, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "kmeans", "stability_kmeans_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "kmeans", "stability_kmeans_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "kmeans", "stability_kmeans_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "kmeans", "stability_kmeans_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




###################
### immunoClust ###
###################

# note installation from Bioconductor requires GNU Scientific Library

library(immunoClust)

source("run_immunoClust_stability.R")


# use non-transformed data files, since immunoClust will transform automatically

file_Levine_32_immunoClust <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransform.fcs")
file_Levine_13_immunoClust <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13_notransform.fcs")
file_Nilsson_immunoClust <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC_notransform.fcs")
file_Mosmann_immunoClust <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ_notransform.fcs")

# input data as flowFrame objects

data_Levine_32_immunoClust <- flowCore::read.FCS(file_Levine_32_immunoClust, transformation = FALSE)
data_Levine_13_immunoClust <- flowCore::read.FCS(file_Levine_13_immunoClust, transformation = FALSE)
data_Nilsson_immunoClust <- flowCore::read.FCS(file_Nilsson_immunoClust, transformation = FALSE)
data_Mosmann_immunoClust <- flowCore::read.FCS(file_Mosmann_immunoClust, transformation = FALSE)

# column names (parameters)

pars_Levine_32 <- colnames(data_Levine_32_immunoClust)[marker_cols_Levine_32]
pars_Levine_13 <- colnames(data_Levine_13_immunoClust)[marker_cols_Levine_13]
pars_Nilsson <- colnames(data_Nilsson_immunoClust)[marker_cols_Nilsson]
pars_Mosmann <- colnames(data_Mosmann_immunoClust)[marker_cols_Mosmann]


# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32_immunoClust, data_Levine_13_immunoClust, 
                           data_Nilsson_immunoClust, data_Mosmann_immunoClust)), n)

res_stability_immunoClust <- bplapply(data_list, run_immunoClust_stability, 
                                      BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_immunoClust, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_immunoClust, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_immunoClust, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_immunoClust, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "immunoClust", "stability_immunoClust_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




#######################
### immunoClust_all ###
#######################

# note installation from Bioconductor requires GNU Scientific Library

library(immunoClust)

source("run_immunoClust_all_stability.R")


# use non-transformed data files, since immunoClust will transform automatically

file_Levine_32_immunoClust_all <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransform.fcs")
file_Levine_13_immunoClust_all <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13_notransform.fcs")
file_Nilsson_immunoClust_all <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC_notransform.fcs")
file_Mosmann_immunoClust_all <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ_notransform.fcs")

# input data as flowFrame objects

data_Levine_32_immunoClust_all <- flowCore::read.FCS(file_Levine_32_immunoClust_all, transformation = FALSE)
data_Levine_13_immunoClust_all <- flowCore::read.FCS(file_Levine_13_immunoClust_all, transformation = FALSE)
data_Nilsson_immunoClust_all <- flowCore::read.FCS(file_Nilsson_immunoClust_all, transformation = FALSE)
data_Mosmann_immunoClust_all <- flowCore::read.FCS(file_Mosmann_immunoClust_all, transformation = FALSE)

# column names (parameters)

pars_Levine_32 <- colnames(data_Levine_32_immunoClust_all)[marker_cols_Levine_32]
pars_Levine_13 <- colnames(data_Levine_13_immunoClust_all)[marker_cols_Levine_13]
pars_Nilsson <- colnames(data_Nilsson_immunoClust_all)[marker_cols_Nilsson]
pars_Mosmann <- colnames(data_Mosmann_immunoClust_all)[marker_cols_Mosmann]


# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32_immunoClust_all, data_Levine_13_immunoClust_all, 
                           data_Nilsson_immunoClust_all, data_Mosmann_immunoClust_all)), n)

res_stability_immunoClust_all <- bplapply(data_list, run_immunoClust_all_stability, 
                                          BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_immunoClust_all, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_immunoClust_all, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_immunoClust_all, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_immunoClust_all, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "immunoClust_all", "stability_immunoClust_all_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")




###################
### SamSPECTRAL ###
###################

library(SamSPECTRAL)

source("run_SamSPECTRAL_stability.R")

# run clustering method n times in parallel

data_list <- rep(list(list(data_Levine_32, data_Levine_13, data_Nilsson, data_Mosmann)), n)

res_stability_SamSPECTRAL <- bplapply(data_list, run_SamSPECTRAL_stability, 
                                      BPPARAM = MulticoreParam(workers = n, RNGseed = 123))

# collapse results into one data frame for each data set

res_stability_Levine_32 <- t(sapply(res_stability_SamSPECTRAL, function(s) s[[1]]))
res_stability_Levine_13 <- t(sapply(res_stability_SamSPECTRAL, function(s) s[[2]]))
res_stability_Nilsson <- t(sapply(res_stability_SamSPECTRAL, function(s) s[[3]]))
res_stability_Mosmann <- t(sapply(res_stability_SamSPECTRAL, function(s) s[[4]]))

# save results to files

write.table(res_stability_Levine_32, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Levine_32.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Levine_13, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Levine_13.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Nilsson, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Nilsson.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_stability_Mosmann, file = file.path(RESULTS_DIR, "SamSPECTRAL", "stability_SamSPECTRAL_Mosmann.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")



