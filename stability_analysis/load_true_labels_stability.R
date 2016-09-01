#########################################################################################
# Stability analysis:
# R script to load true cluster labels (manually gated populations)
#
# Lukas Weber, August 2016
#########################################################################################


library(flowCore)

DATA_DIR <- "../../../benchmark_data_sets"


# load true cluster labels (manually gated populations)

files_truth <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare.fcs")
)

clus_truth <- vector("list", length(files_truth))
names(clus_truth) <- names(files_truth)

for (i in 1:length(clus_truth)) {
  data_truth_i <- flowCore::exprs(flowCore::read.FCS(files_truth[[i]], transformation = FALSE, truncate_max_range = FALSE))
  clus_truth[[i]] <- data_truth_i[, "label"]
}

sapply(clus_truth, length)
sapply(clus_truth, table)

