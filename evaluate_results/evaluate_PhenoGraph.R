#########################################################################################
# R script to load and evaluate results for PhenoGraph
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
RES_DIR_PHENOGRAPH <- "../../results_auto/PhenoGraph"

DATA_DIR <- "../../../benchmark_data_sets"

# which data sets required subsampling for this method (see parameters spreadsheet)
is_subsampled <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

is_rare <- c(FALSE, FALSE, FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE)

# note: no FlowCAP data sets for this method




####################################################
### load truth (manual gating population labels) ###
####################################################

# files with true population labels (subsampled labels if subsampling was required for
# this method; see parameters spreadsheet)

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
  if (!is_subsampled[i]) {
    data_truth_i <- flowCore::exprs(flowCore::read.FCS(files_truth[[i]], transformation = FALSE, truncate_max_range = FALSE))
  } else {
    data_truth_i <- read.table(files_truth[[i]], header = TRUE, stringsAsFactors = FALSE)
  }
  clus_truth[[i]] <- data_truth_i[, "label"]
}

sapply(clus_truth, length)

# cluster sizes and number of clusters
# (for data sets with single rare population: 1 = rare population of interest, 0 = all others)

tbl_truth <- lapply(clus_truth, table)

tbl_truth
sapply(tbl_truth, length)




###############################
### load PhenoGraph results ###
###############################

# load cluster labels

files_out <- list(
  Levine_32dim = file.path(RES_DIR_PHENOGRAPH, "PhenoGraph_results_Levine_32dim.fcs"), 
  Levine_13dim = file.path(RES_DIR_PHENOGRAPH, "PhenoGraph_results_Levine_13dim.fcs"), 
  Samusik_01   = file.path(RES_DIR_PHENOGRAPH, "PhenoGraph_results_Samusik_01.fcs"), 
  Samusik_all  = file.path(RES_DIR_PHENOGRAPH, "PhenoGraph_results_Samusik_all.fcs"), 
  Nilsson_rare = file.path(RES_DIR_PHENOGRAPH, "PhenoGraph_results_Nilsson_rare.fcs"), 
  Mosmann_rare = file.path(RES_DIR_PHENOGRAPH, "PhenoGraph_results_Mosmann_rare.fcs")
)

clus <- vector("list", length(files_out))
names(clus) <- names(files_out)

for (i in 1:length(clus)) {
  data_i <- flowCore::exprs(flowCore::read.FCS(files_out[[i]], transformation = FALSE, truncate_max_range = FALSE))
  clus[[i]] <- data_i[, grep("PhenoGraph", colnames(data_i))]
}

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

# return named object (used in plotting scripts)

res_PhenoGraph <- res



