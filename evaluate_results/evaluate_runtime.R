#########################################################################################
# R script to load runtimes
#
# Lukas Weber, September 2016
#########################################################################################

# Methods that were run using R scripts were timed with the "system.time()" function, 
# with results saved as .txt files in folders "results/auto/runtime/" or 
# "results/manual/runtime/".

# Methods not accessible via R were timed manually using a stopwatch app or the Matlab 
# "clock" command, with results saved in the spreadsheet 
# "Supp_Table_S1_clustering_methods_parameters.xlsx". These values need to be entered
# manually here.



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
# Manually enter runtimes for methods not accessible via R (values from spreadsheet
# "clustering_methods_details.xlsx")
# =================================================================================

res_runtime[["Levine_32dim"]][["ACCENSE"]] <- (5 * 60) + 32
res_runtime[["Levine_13dim"]][["ACCENSE"]] <- (4 * 60) + 48
res_runtime[["Samusik_01"]][["ACCENSE"]]   <- (6 * 60) + 21
res_runtime[["Samusik_all"]][["ACCENSE"]]  <- (5 * 60) + 32
res_runtime[["Nilsson_rare"]][["ACCENSE"]] <- (6 * 60) + 11
res_runtime[["Mosmann_rare"]][["ACCENSE"]] <- (4 * 60) + 37

res_runtime[["Levine_32dim"]][["PhenoGraph"]] <- (37 * 60) + 0
res_runtime[["Levine_13dim"]][["PhenoGraph"]] <- (12 * 60) + 09
res_runtime[["Samusik_01"]][["PhenoGraph"]]   <- (5 * 60) + 55
res_runtime[["Samusik_all"]][["PhenoGraph"]]  <- (5 * 3600) + (30 * 60) + 35
res_runtime[["Nilsson_rare"]][["PhenoGraph"]] <- (1 * 60) + 58
res_runtime[["Mosmann_rare"]][["PhenoGraph"]] <- (43 * 60) + 43

res_runtime[["Levine_32dim"]][["SWIFT"]] <- (2 * 3600) + (27 * 60) + 39
res_runtime[["Levine_13dim"]][["SWIFT"]] <- (1 * 3600) + (7 * 60) + 3
res_runtime[["Samusik_01"]][["SWIFT"]]   <- (2 * 3600) + (19 * 60) + 30
res_runtime[["Samusik_all"]][["SWIFT"]]  <- (2 * 3600) + (50 * 60) + 8
# note using results from SWIFT version 2.0 for Nilsson_rare and Mosmann_rare
res_runtime[["Nilsson_rare"]][["SWIFT"]] <- (11 * 60) + 26
res_runtime[["Mosmann_rare"]][["SWIFT"]] <- (34 * 60) + 34

res_runtime[["Levine_32dim"]][["Xshift"]] <- (4 * 3600) + (45 * 60) + 26
res_runtime[["Levine_13dim"]][["Xshift"]] <- (48 * 60) + 17
res_runtime[["Samusik_01"]][["Xshift"]]   <- (24 * 60) + 54
res_runtime[["Samusik_all"]][["Xshift"]]  <- (3 * 3600) + (48 * 60) + 27
res_runtime[["Nilsson_rare"]][["Xshift"]] <- (4 * 60) + 37
res_runtime[["Mosmann_rare"]][["Xshift"]] <- (3 * 3600) + (18 * 60) + 20



# ===================================
# Load runtimes for all other methods
# ===================================

# run script "evaluate_all_methods.R" to load results directories (automatic or manual
# number of clusters)

#source("evaluate_all_methods.R")  ## takes 20 min

RES_DIRS <- c(RES_DIR_CLUSTERX, 
              RES_DIR_DENSVM, 
              RES_DIR_FLOCK, 
              RES_DIR_FLOWCLUST, 
              RES_DIR_FLOWMEANS, 
              RES_DIR_FLOWMERGE, 
              RES_DIR_FLOWPEAKS, 
              RES_DIR_FLOWSOM, 
              RES_DIR_FLOWSOM_PRE, 
              RES_DIR_IMMUNOCLUST, 
              RES_DIR_KMEANS, 
              RES_DIR_RCLUSTERPP, 
              RES_DIR_SAMSPECTRAL, 
              RES_DIR_SPADE)

method_names <- c("ClusterX", 
                  "DensVM", 
                  "FLOCK", 
                  "flowClust", 
                  "flowMeans", 
                  "flowMerge", 
                  "flowPeaks", 
                  "FlowSOM", 
                  "FlowSOM_pre", 
                  "immunoClust", 
                  "kmeans", 
                  "Rclusterpp", 
                  "SamSPECTRAL", 
                  "SPADE")


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



# ===========
# Subsampling
# ===========

# which methods required subsampling (see parameters spreadsheet); not including methods 
# that were not included in final results (e.g. due to errors or non-completion)

which_sub_Levine_32dim <- c("ACCENSE", "ClusterX", "DensVM", "immunoClust", "SWIFT")  ## exclude: flowClust, flowMerge
which_sub_Levine_13dim <- c("ACCENSE", "DensVM", "flowClust", "flowMerge")
which_sub_Samusik_01   <- c("ACCENSE", "flowClust", "flowMerge")
which_sub_Samusik_all  <- c("ACCENSE", "ClusterX", "DensVM", "flowClust", "flowMeans", "flowMerge", 
                            "immunoClust", "Rclusterpp", "SamSPECTRAL", "SWIFT", "Xshift")
which_sub_Nilsson_rare <- c("ACCENSE", "flowMerge")
which_sub_Mosmann_rare <- c("ACCENSE", "ClusterX", "DensVM", "flowClust", "flowMerge", "Rclusterpp")
which_sub_FlowCAP_ND   <- c("ClusterX", "DensVM", "flowMerge")
which_sub_FlowCAP_WNV  <- c("ClusterX", "DensVM")



# ========================
# Multiple processor cores
# ========================

# which methods required multiple processor cores (see parameters spreadsheeet); not
# including methods that were not included in final results (e.g. due to errors)

which_cores_Levine_32dim <- c("Rclusterpp", "SWIFT", "Xshift")  ## exclude: SPADE
which_cores_Levine_13dim <- c("Rclusterpp", "SPADE", "SWIFT", "Xshift")
which_cores_Samusik_01   <- c("Rclusterpp", "SPADE", "SWIFT", "Xshift")
which_cores_Samusik_all  <- c("Rclusterpp", "SPADE", "SWIFT", "Xshift")
which_cores_Nilsson_rare <- c("Rclusterpp", "SPADE", "SWIFT", "Xshift")
which_cores_Mosmann_rare <- c("Rclusterpp", "SPADE", "SWIFT", "Xshift")
which_cores_FlowCAP_ND   <- c("Rclusterpp", "SPADE")  ## exclude: SWIFT, X-shift
which_cores_FlowCAP_WNV  <- c("Rclusterpp")  ## exclude: SPADE, SWIFT, X-shift


