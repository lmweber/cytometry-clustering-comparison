#########################################################################################
# R script to load runtime results
#
# Lukas M. Weber, March 2016
#########################################################################################


# Methods that were not available as R packages were timed manually with a stopwatch app, 
# with results saved in the spreadsheet "clustering_methods_details.xlsx".

# All other methods were timed with system.time() in R, with results saved as text files
# in the folder "results_auto_n_clus/runtime/" or "results_manual_n_clus/runtime/".



# ===================
# Prepare data frames
# ===================

runtime_Levine_32 <- list()
runtime_Levine_13 <- list()
runtime_Nilsson <- list()
runtime_Mosmann <- list()



# ===================================================================================
# Manually enter runtime results for methods not available as R packages (copied from
# spreadsheet "clustering_methods_details.xlsx")
# ===================================================================================

runtime_Levine_32[["ACCENSE"]] <- (4 * 60) + 49
runtime_Levine_13[["ACCENSE"]] <- (4 * 60) + 41
runtime_Nilsson[["ACCENSE"]] <- (5 * 60) + 40
runtime_Mosmann[["ACCENSE"]] <- (4 * 60) + 7

runtime_Levine_32[["FLOCK"]] <- (2 * 60) + 51
runtime_Levine_13[["FLOCK"]] <- 24
runtime_Nilsson[["FLOCK"]] <- 6
runtime_Mosmann[["FLOCK"]] <- (1 * 60) + 1

runtime_Levine_32[["PhenoGraph"]] <- (37 * 60) + 42
runtime_Levine_13[["PhenoGraph"]] <- (13 * 60) + 11
runtime_Nilsson[["PhenoGraph"]] <- (2 * 60) + 2
runtime_Mosmann[["PhenoGraph"]] <- (49 * 60) + 59

runtime_Levine_32[["SWIFT"]] <- (2 * 3600) + (35 * 60) + 48
runtime_Levine_13[["SWIFT"]] <- (56 * 60) + 13
runtime_Nilsson[["SWIFT"]] <- (11 * 60) + 26
runtime_Mosmann[["SWIFT"]] <- (34 * 60) + 34



# ==========================================
# Load runtime results for all other methods
# ==========================================

# load results directories (depending on whether automatic or manually selected number of clusters)

source("load_results_directories.R")

RES_DIRS <- c(RES_DIR_DENSVM, 
              RES_DIR_FLOWMEANS, 
              RES_DIR_FLOWSOM, 
              RES_DIR_FLOWSOM_META, 
              RES_DIR_IMMUNOCLUST, 
              RES_DIR_IMMUNOCLUST_ALL, 
              RES_DIR_KMEANS, 
              RES_DIR_RCLUSTERPP, 
              RES_DIR_SAMSPECTRAL)

method_names <- c("DensVM", 
                  "flowMeans", 
                  "FlowSOM", 
                  "FlowSOM_meta", 
                  "immunoClust", 
                  "immunoClust_all", 
                  "kmeans", 
                  "Rclusterpp", 
                  "SamSPECTRAL")


# load runtime results files

for (i in 1:length(RES_DIRS)) {
  file_i <- paste0(RES_DIRS[i], "/runtime/runtime_", method_names[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  runtime_Levine_32[[method_names[i]]] <- data_i["Levine_2015_marrow_32", "runtime"]
  runtime_Levine_13[[method_names[i]]] <- data_i["Levine_2015_marrow_13", "runtime"]
  runtime_Nilsson[[method_names[i]]] <- data_i["Nilsson_2013_HSC", "runtime"]
  runtime_Mosmann[[method_names[i]]] <- data_i["Mosmann_2014_activ", "runtime"]
}



# ==========================
# Arrange in ascending order
# ==========================

runtime_Levine_32_ord <- unlist(runtime_Levine_32)
runtime_Levine_13_ord <- unlist(runtime_Levine_13)
runtime_Nilsson_ord <- unlist(runtime_Nilsson)
runtime_Mosmann_ord <- unlist(runtime_Mosmann)

runtime_Levine_32_ord <- runtime_Levine_32_ord[order(runtime_Levine_32_ord)]
runtime_Levine_13_ord <- runtime_Levine_13_ord[order(runtime_Levine_13_ord)]
runtime_Nilsson_ord <- runtime_Nilsson_ord[order(runtime_Nilsson_ord)]
runtime_Mosmann_ord <- runtime_Mosmann_ord[order(runtime_Mosmann_ord)]


