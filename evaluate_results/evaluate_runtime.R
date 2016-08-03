#########################################################################################
# R script to load runtimes
#
# Lukas Weber, August 2016
#########################################################################################

# Methods that were run using R scripts were timed with the "system.time()" function, 
# with results saved as .txt files in folders "results_auto/runtime/" or 
# "results_manual/runtime/".

# Methods not accessible via R were timed manually using a stopwatch app or the Matlab
# "clock" command, with results saved in the spreadsheet 
# "Supp_Table_S1_clustering_methods_parameters.xlsx".



# ===============
# Prepare objects
# ===============

res_runtime <- list(
  Levine_32dim = list(), 
  Levine_13dim = list(), 
  Samusik_01   = list(), 
  Samusik_all  = list(), 
  Nilsson_rare = list(), 
  Mosmann_rare = list(), 
  FlowCAP_ND   = list(), 
  FlowCAP_WNV  = list()
)



# =================================================================================
# Manually enter runtimes for methods not accessible via R (copied from spreadsheet
# "clustering_methods_details.xlsx")
# =================================================================================

res_runtime[["Levine_32dim"]][["ACCENSE"]] <- (4 * 60) + 49
res_runtime[["Levine_13dim"]][["ACCENSE"]] <- (4 * 60) + 41
res_runtime[["Nilsson_rare"]][["ACCENSE"]] <- (5 * 60) + 40
res_runtime[["Mosmann_rare"]][["ACCENSE"]] <- (4 * 60) + 7

res_runtime[["Levine_32dim"]][["PhenoGraph"]] <- (37 * 60) + 42
res_runtime[["Levine_13dim"]][["PhenoGraph"]] <- (13 * 60) + 11
res_runtime[["Nilsson_rare"]][["PhenoGraph"]] <- (2 * 60) + 2
res_runtime[["Mosmann_rare"]][["PhenoGraph"]] <- (49 * 60) + 59

res_runtime[["Levine_32dim"]][["SWIFT"]] <- (2 * 3600) + (35 * 60) + 48
res_runtime[["Levine_13dim"]][["SWIFT"]] <- (56 * 60) + 13
res_runtime[["Nilsson_rare"]][["SWIFT"]] <- (11 * 60) + 26
res_runtime[["Mosmann_rare"]][["SWIFT"]] <- (34 * 60) + 34



# ===================================
# Load runtimes for all other methods
# ===================================

# run script "evaluate_all_methods.R" to load results directories (automatic or manual
# number of clusters)

#source("evaluate_all_methods.R")  ## takes 15 min

RES_DIRS <- c(RES_DIR_CLUSTERX, 
              RES_DIR_DENSVM, 
              RES_DIR_FLOCK, 
              RES_DIR_FLOWMEANS, 
              RES_DIR_FLOWPEAKS, 
              RES_DIR_FLOWSOM, 
              RES_DIR_FLOWSOM_PRE, 
              RES_DIR_IMMUNOCLUST, 
              RES_DIR_KMEANS, 
              RES_DIR_RCLUSTERPP)

method_names <- c("ClusterX", 
                  "DensVM", 
                  "FLOCK", 
                  "flowMeans", 
                  "flowPeaks", 
                  "FlowSOM", 
                  "FlowSOM_pre", 
                  "immunoClust", 
                  "kmeans", 
                  "Rclusterpp")


# load runtime results files

for (i in 1:length(RES_DIRS)) {
  file_i <- paste0(RES_DIRS[i], "/../runtimes/runtime_", method_names[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  res_runtime[["Levine_32dim"]][[method_names[i]]] <- data_i["Levine_32dim", "runtime"]
  res_runtime[["Levine_13dim"]][[method_names[i]]] <- data_i["Levine_13dim", "runtime"]
  res_runtime[["Samusik_01"]][[method_names[i]]]   <- data_i["Samusik_01", "runtime"]
  res_runtime[["Samusik_all"]][[method_names[i]]]  <- data_i["Samusik_all", "runtime"]
  res_runtime[["Nilsson_rare"]][[method_names[i]]] <- data_i["Nilsson_rare", "runtime"]
  res_runtime[["Mosmann_rare"]][[method_names[i]]] <- data_i["Mosmann_rare", "runtime"]
  res_runtime[["FlowCAP_ND"]][[method_names[i]]]   <- data_i["FlowCAP_ND", "runtime"]
  res_runtime[["FlowCAP_WNV"]][[method_names[i]]]  <- data_i["FlowCAP_WNV", "runtime"]
}



# ==========================
# Arrange in ascending order
# ==========================

res_runtime_ord <- lapply(res_runtime, function(r) {
  r_ord <- unlist(r)
  r_ord[order(r_ord)]
})


