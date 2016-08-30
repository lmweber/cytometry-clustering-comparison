#########################################################################################
# Stability analysis (multiple random starts):
# R script to run and evaluate all methods. Methods are run multiple times in parallel 
# using BiocParallel::bplapply() (where possible).
# 
# Lukas Weber, August 2016
#########################################################################################


library(flowCore)
library(clue)
library(BiocParallel)

library(flowMeans)
library(flowPeaks)
library(FlowSOM)
library(immunoClust)
library(Rclusterpp)
library(SamSPECTRAL)

# helper functions to match clusters and evaluate
source("../helpers/helper_match_evaluate_multiple.R")
source("../helpers/helper_match_evaluate_single.R")

# load data and true population labels
source("stability_load_data.R")
source("stability_load_true_labels.R")

# functions to run and evaluate methods
source("random_starts_FLOCK.R")
source("random_starts_flowMeans.R")
source("random_starts_flowPeaks.R")
source("random_starts_FlowSOM_pre.R")
source("random_starts_FlowSOM.R")
source("random_starts_immunoClust.R")
source("random_starts_kmeans.R")
source("random_starts_Rclusterpp.R")
source("random_starts_SamSPECTRAL.R")

# directory to save results
RESULTS_DIR <- "../../results_stability_random_starts"

# number of times to run each method
n <- 30

# replicate data
data <- rep(list(data), n)
data_notransform <- rep(list(data_notransform), n)




###########################################
### RUN AND EVALUATE METHODS (PARALLEL) ###
###########################################

# run methods and evaluate; using BiocParallel::bplapply() for parallelization

seed <- 123

res_random_starts_flowMeans <- bplapply(data, random_starts_flowMeans, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : flowMeans complete\n")

res_random_starts_flowPeaks <- bplapply(data, random_starts_flowPeaks, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : flowPeaks complete\n")

res_random_starts_FlowSOM_pre <- bplapply(data, random_starts_FlowSOM_pre, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : FlowSOM_pre complete\n")

res_random_starts_FlowSOM <- bplapply(data, random_starts_FlowSOM, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : FlowSOM complete\n")

# immunoClust: data set Mosmann_rare only (due to subsampling); also use non-transformed data due to automatic transform
data_immunoClust <- lapply(data_notransform, function(l) l["Mosmann_rare"])
res_random_starts_immunoClust <- bplapply(data_immunoClust, random_starts_immunoClust, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : immunoClust complete\n")

res_random_starts_kmeans <- bplapply(data, random_starts_kmeans, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : kmeans complete\n")

# Rclusterpp: data set Levine_32dim only (due to subsampling)
data_Rclusterpp <- lapply(data, function(l) l["Levine_32dim"])
res_random_starts_Rclusterpp <- bplapply(data_Rclusterpp, random_starts_Rclusterpp, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : Rclusterpp complete\n")

res_random_starts_SamSPECTRAL <- bplapply(data, random_starts_SamSPECTRAL, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : SamSPECTRAL complete\n")




#########################################
### RUN AND EVALUATE METHODS (SERIES) ###
#########################################

# run methods multiple times in series (parallelization not possible)

set.seed(seed)
res_random_starts_FLOCK <- lapply(data, random_starts_FLOCK)
cat("stability analysis (random starts) : FLOCK complete\n")




####################
### SAVE RESULTS ###
####################

res_random_starts <- list(FLOCK = res_random_starts_FLOCK, 
                          flowMeans = res_random_starts_flowMeans, 
                          flowPeaks = res_random_starts_flowPeaks, 
                          FlowSOM_pre = res_random_starts_FlowSOM_pre, 
                          FlowSOM = res_random_starts_FlowSOM, 
                          immunoClust = res_random_starts_immunoClust, 
                          kmeans = res_random_starts_kmeans, 
                          Rclusterpp = res_random_starts_Rclusterpp, 
                          SamSPECTRAL = res_random_starts_SamSPECTRAL)

# remove any methods skipped for each data set
res_Levine_32dim <- res_random_starts[-which(names(res_random_starts) == "immunoClust")]
res_Mosmann_rare <- res_random_starts[-which(names(res_random_starts) == "Rclusterpp")]


# collapse into one data frame per data set

res_Levine_32dim <- lapply(res_Levine_32dim, 
                           function(r) t(sapply(r, function(l) l[["Levine_32dim"]][c("mean_pr", "mean_re", "mean_F1")])))

res_Mosmann_rare <- lapply(res_Mosmann_rare, 
                           function(r) t(sapply(r, function(l) l[["Mosmann_rare"]][c("pr", "re", "F1")])))


# save files

for (i in 1:length(res_Levine_32dim)) {
  method_name <- names(res_Levine_32dim)[i]
  file = file.path(RESULTS_DIR, method_name, paste0("stability_random_starts_", method_name, "_Levine_32dim.txt"))
  write.table(res_Levine_32dim[[i]], file = file, row.names = FALSE, quote = FALSE, sep = "\t")
}

for (i in 1:length(res_Mosmann_rare)) {
  method_name <- names(res_Mosmann_rare)[i]
  file = file.path(RESULTS_DIR, method_name, paste0("stability_random_starts_", method_name, "_Mosmann_rare.txt"))
  write.table(res_Mosmann_rare[[i]], file = file, row.names = FALSE, quote = FALSE, sep = "\t")
}



