#########################################################################################
# R script to load and evaluate results for ACCENSE
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
RES_DIR_ACCENSE <- "../../results/auto/ACCENSE"

DATA_DIR <- "../../../benchmark_data_sets"

# which data sets required subsampling for this method (see parameters spreadsheet)
is_subsampled <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

is_rare <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)

# note: no FlowCAP data sets for this method




####################################################
### load truth (manual gating population labels) ###
####################################################

# files with true population labels (subsampled labels if subsampling was required for
# this method; see parameters spreadsheet)

files_truth <- list(
  Levine_32dim = file.path(RES_DIR_ACCENSE, "accense_output_Levine_32dim.csv"), 
  Levine_13dim = file.path(RES_DIR_ACCENSE, "accense_output_Levine_13dim.csv"), 
  Samusik_01   = file.path(RES_DIR_ACCENSE, "accense_output_Samusik_01.csv"), 
  Samusik_all  = file.path(RES_DIR_ACCENSE, "accense_output_Samusik_all.csv"), 
  Nilsson_rare = file.path(RES_DIR_ACCENSE, "accense_output_Nilsson_rare.csv"), 
  Mosmann_rare = file.path(RES_DIR_ACCENSE, "accense_output_Mosmann_rare.csv")
)

# extract true population labels

clus_truth <- vector("list", length(files_truth))
names(clus_truth) <- names(files_truth)

for (i in 1:length(clus_truth)) {
  if (!is_subsampled[i]) {
    data_truth_i <- flowCore::exprs(flowCore::read.FCS(files_truth[[i]], transformation = FALSE, truncate_max_range = FALSE))
  } else if (is_subsampled[i]) {
    data_truth_i <- read.csv(files_truth[[i]], stringsAsFactors = FALSE)
  }
  clus_truth[[i]] <- data_truth_i[, "label"]
}

sapply(clus_truth, length)

# cluster sizes and number of clusters
# (for data sets with single rare population: 1 = rare population of interest, 0 = all others)

tbl_truth <- lapply(clus_truth, table)

tbl_truth
sapply(tbl_truth, length)

# store named objects (for other scripts)

files_truth_ACCENSE <- files_truth
clus_truth_ACCENSE <- clus_truth




############################
### load ACCENSE results ###
############################

# load cluster labels

files_out <- files_truth

clus <- lapply(files_out, function(f) {
  read.csv(f, stringsAsFactors = FALSE)[, "population"]
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

files_ACCENSE <- files_out
clus_ACCENSE <- clus




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

res_ACCENSE <- res



