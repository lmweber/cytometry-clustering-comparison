#########################################################################################
# R script to load and evaluate results for ensemble clustering
#
# Lukas Weber, September 2016
#########################################################################################


library(flowCore)
library(clue)

# helper functions to match clusters and evaluate
source("../helpers/helper_match_evaluate_multiple.R")
source("../helpers/helper_match_evaluate_single.R")

# directory with saved results
RES_DIR_ENSEMBLE <- "../../results/ensemble"

DATA_DIR <- "../../../benchmark_data_sets"

# alternate FlowCAP results at the end
is_rare    <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)




####################################################
### load truth (manual gating population labels) ###
####################################################

# files with true population labels

files_truth <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik/data/Samusik_all.fcs"), 
  Nilsson_rare = file.path(DATA_DIR, "Nilsson_rare/data/Nilsson_rare.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare.fcs")
)

# extract true population labels

clus_truth <- vector("list", length(files_truth))
names(clus_truth) <- names(files_truth)

for (i in 1:length(clus_truth)) {
  data_truth_i <- flowCore::exprs(flowCore::read.FCS(files_truth[[i]], transformation = FALSE, truncate_max_range = FALSE))
  clus_truth[[i]] <- data_truth_i[, "label"]
}

sapply(clus_truth, length)

# cluster sizes and number of clusters
# (for data sets with single rare population: 1 = rare population of interest, 0 = all others)

tbl_truth <- lapply(clus_truth, table)

tbl_truth
sapply(tbl_truth, length)




########################################
### load ensemble clustering results ###
########################################

# load cluster labels

files_out <- list(
  Levine_32dim = file.path(RES_DIR_ENSEMBLE, "ensemble_labels_Levine_32dim.txt"), 
  Levine_13dim = file.path(RES_DIR_ENSEMBLE, "ensemble_labels_Levine_13dim.txt"), 
  Samusik_01   = file.path(RES_DIR_ENSEMBLE, "ensemble_labels_Samusik_01.txt"), 
  Samusik_all  = file.path(RES_DIR_ENSEMBLE, "ensemble_labels_Samusik_all.txt"), 
  Nilsson_rare = file.path(RES_DIR_ENSEMBLE, "ensemble_labels_Nilsson_rare.txt"), 
  Mosmann_rare = file.path(RES_DIR_ENSEMBLE, "ensemble_labels_Mosmann_rare.txt")
)

clus <- lapply(files_out, function(f) {
  read.table(f, header = TRUE, stringsAsFactors = FALSE)[, "label"]
})

sapply(clus, length)

# cluster sizes and number of clusters
# (for data sets with single rare population: 1 = rare population of interest, 0 = all others)

tbl <- lapply(clus, table)

tbl
sapply(tbl, length)

# contingency tables

for (i in 1:length(clus)) {
  print(table(clus[[i]], clus_truth[[i]]))
}

# store named object

clus_ensemble <- clus




###################################
### match clusters and evaluate ###
###################################

# see helper function scripts for details on matching strategy and evaluation

res <- vector("list", length(clus))
names(res) <- names(clus)

for (i in 1:length(clus)) {
  if (!is_rare[i]) {
    res[[i]] <- helper_match_evaluate_multiple(clus[[i]], clus_truth[[i]])
    
  } else if (is_rare[i]) {
    res[[i]] <- helper_match_evaluate_single(clus[[i]], clus_truth[[i]])
  }
}

# store named object (for plotting scripts)

res_ensemble <- res



