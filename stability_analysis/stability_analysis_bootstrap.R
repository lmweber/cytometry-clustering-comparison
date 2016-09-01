#########################################################################################
# Stability analysis (bootstrap resamples):
# R script to run and evaluate all methods. Methods are run multiple times in parallel 
# using BiocParallel::bplapply() (where possible).
# 
# Lukas Weber, September 2016
#########################################################################################


library(flowCore)
library(clue)
library(BiocParallel)

library(flowMeans)
library(flowPeaks)
library(FlowSOM)
library(immunoClust)
library(SamSPECTRAL)

# helper functions to match clusters and evaluate
source("../helpers/helper_match_evaluate_multiple.R")
source("../helpers/helper_match_evaluate_single.R")

# load data and true population labels
source("load_data_stability.R")
source("load_true_labels_stability.R")

# functions to run and evaluate methods
source("run_FLOCK_stability.R")
source("run_flowMeans_stability.R")
source("run_flowPeaks_stability.R")
source("run_FlowSOM_pre_stability.R")
source("run_FlowSOM_stability.R")
source("run_immunoClust_stability.R")
source("run_kmeans_stability.R")
source("run_SamSPECTRAL_stability.R")

# directory to save results
RESULTS_DIR <- "../../results_stability_bootstrap"

# number of times to run each method
n <- 30




###########################################
### REPLICATE DATA: BOOTSTRAP RESAMPLES ###
###########################################

# replicate data n times (bootstrap resamples)

fn_bootstrap <- function(d) {
  d <- d[sample(1:nrow(d), nrow(d), replace = TRUE), ]
}

set.seed(123)
data <- lapply(1:n, function(i) lapply(data, fn_bootstrap))

set.seed(123)
data_notransform <- lapply(1:n, function(i) lapply(data_notransform, fn_bootstrap))




###########################################
### RUN AND EVALUATE METHODS (PARALLEL) ###
###########################################

# run methods and evaluate; using BiocParallel::bplapply() for parallelization

seed <- 123

res_bootstrap_flowMeans <- bplapply(data, run_flowMeans_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : flowMeans complete\n")

res_bootstrap_flowPeaks <- bplapply(data, run_flowPeaks_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : flowPeaks complete\n")

res_bootstrap_FlowSOM_pre <- bplapply(data, run_FlowSOM_pre_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : FlowSOM_pre complete\n")

res_bootstrap_FlowSOM <- bplapply(data, run_FlowSOM_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : FlowSOM complete\n")

# immunoClust: data set Mosmann_rare only (due to subsampling); also use non-transformed data due to automatic transform
data_immunoClust <- lapply(data_notransform, function(l) l["Mosmann_rare"])
res_bootstrap_immunoClust <- bplapply(data_immunoClust, run_immunoClust_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : immunoClust complete\n")

res_bootstrap_kmeans <- bplapply(data, run_kmeans_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : kmeans complete\n")

res_bootstrap_SamSPECTRAL <- bplapply(data, run_SamSPECTRAL_stability, BPPARAM = MulticoreParam(workers = n, RNGseed = seed))
cat("stability analysis (random starts) : SamSPECTRAL complete\n")




#########################################
### RUN AND EVALUATE METHODS (SERIES) ###
#########################################

# run methods multiple times in series (parallelization not possible)

set.seed(seed)
res_bootstrap_FLOCK <- lapply(data, run_FLOCK_stability)
cat("stability analysis (random starts) : FLOCK complete\n")




####################
### SAVE RESULTS ###
####################

res_bootstrap <- list(FLOCK = res_bootstrap_FLOCK, 
                      flowMeans = res_bootstrap_flowMeans, 
                      flowPeaks = res_bootstrap_flowPeaks, 
                      FlowSOM_pre = res_bootstrap_FlowSOM_pre, 
                      FlowSOM = res_bootstrap_FlowSOM, 
                      immunoClust = res_bootstrap_immunoClust, 
                      kmeans = res_bootstrap_kmeans, 
                      SamSPECTRAL = res_bootstrap_SamSPECTRAL)

# remove any methods skipped for each data set
res_Levine_32dim <- res_bootstrap[-which(names(res_bootstrap) == "immunoClust")]
res_Mosmann_rare <- res_bootstrap


# collapse into one data frame per data set

res_Levine_32dim <- lapply(res_Levine_32dim, 
                           function(r) t(sapply(r, function(l) l[["Levine_32dim"]][c("mean_pr", "mean_re", "mean_F1")])))

res_Mosmann_rare <- lapply(res_Mosmann_rare, 
                           function(r) t(sapply(r, function(l) l[["Mosmann_rare"]][c("pr", "re", "F1")])))


# save files

for (i in 1:length(res_Levine_32dim)) {
  method_name <- names(res_Levine_32dim)[i]
  file = file.path(RESULTS_DIR, method_name, paste0("stability_bootstrap_", method_name, "_Levine_32dim.txt"))
  write.table(res_Levine_32dim[[i]], file = file, row.names = FALSE, quote = FALSE, sep = "\t")
}

for (i in 1:length(res_Mosmann_rare)) {
  method_name <- names(res_Mosmann_rare)[i]
  file = file.path(RESULTS_DIR, method_name, paste0("stability_bootstrap_", method_name, "_Mosmann_rare.txt"))
  write.table(res_Mosmann_rare[[i]], file = file, row.names = FALSE, quote = FALSE, sep = "\t")
}



