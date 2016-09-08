#########################################################################################
# R script to load and evaluate results for SWIFT
#
# Lukas Weber, August 2016
#########################################################################################


library(flowCore)
library(clue)

# helper functions to match clusters and evaluate
source("../helpers/helper_match_evaluate_multiple.R")
source("../helpers/helper_match_evaluate_single.R")
source("../helpers/helper_match_evaluate_FlowCAP.R")
source("../helpers/helper_match_evaluate_FlowCAP_alternate.R")

# which set of results to use: automatic or manual number of clusters (see parameters spreadsheet)
RES_DIR_SWIFT <- "../../results/auto/SWIFT"

DATA_DIR <- "../../../benchmark_data_sets"

# which data sets required subsampling for this method (see parameters spreadsheet)
is_subsampled <- c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

is_rare <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)

# note: no FlowCAP data sets for this method




####################################################
### load truth (manual gating population labels) ###
####################################################

# files with true population labels (subsampled labels if subsampling was required for
# this method; see parameters spreadsheet)

files_truth <- list(
  Levine_32dim = file.path(RES_DIR_SWIFT, "Levine_32dim_notransform_subsampled.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01.fcs"), 
  Samusik_all  = file.path(RES_DIR_SWIFT, "Samusik_all_notransform_subsampled.fcs"), 
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

# store named objects (for other scripts)

files_truth_SWIFT <- files_truth
clus_truth_SWIFT <- clus_truth




##########################
### load SWIFT results ###
##########################

# note: using results from older version of SWIFT (2.0) for data sets Nilsson_rare and 
# Mosmann_rare -- the latest version (3.0) has an error during data filtering, which 
# results in almost all cells being removed from these data sets (and data filtering 
# cannot be disabled)

# load cluster labels

files_out <- list(
  Levine_32dim = file.path(RES_DIR_SWIFT, "Levine_32dim_notransform_subsampled.fcs.Cluster_Output.txt"), 
  Levine_13dim = file.path(RES_DIR_SWIFT, "Levine_13dim_notransform.fcs.Cluster_Output.txt"), 
  Samusik_01   = file.path(RES_DIR_SWIFT, "Samusik_01_notransform.fcs.Cluster_Output.txt"), 
  Samusik_all  = file.path(RES_DIR_SWIFT, "Samusik_all_notransform_subsampled.fcs.Cluster_Output.txt"), 
  Nilsson_rare = file.path(RES_DIR_SWIFT, "SWIFT_version_2.0/Nilsson_2013_HSC/Nilsson_2013_HSC_notransform.fcs.Cluster_Output.txt"), 
  Mosmann_rare = file.path(RES_DIR_SWIFT, "SWIFT_version_2.0/Mosmann_2014_activ/Mosmann_2014_activ_notransform.fcs.Cluster_Output.txt")
)

clus <- lapply(files_out, function(f) {
  read.table(f, header = TRUE, sep = "\t", comment.char = "")[, "MergeCluster."]
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

# store named objects (for other scripts)

files_SWIFT <- files_out
clus_SWIFT <- clus




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

res_SWIFT <- res



