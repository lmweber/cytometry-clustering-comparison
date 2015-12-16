#########################################################################################
# R script to load runtime results
#
# Lukas M. Weber, December 2015
#########################################################################################


# Methods that were not available as R packages were timed manually with a stopwatch, 
# with results saved in the spreadsheet "parameters_and_runtime.xlsx".

# All other methods were timed with system.time() in R, with results saved as text files
# in the folder "results/runtime".



# ===================
# Prepare data frames
# ===================

runtime_Levine_32 <- list()
runtime_Levine_13 <- list()
runtime_Nilsson <- list()
runtime_Mosmann <- list()



# ===================================================================================
# Manually enter runtime results for methods not available as R packages (copied from
# spreadsheet "parameters_and_runtime.xlsx")
# ===================================================================================

runtime_Levine_32[["ACCENSE"]] <- (5 * 60) + 12
runtime_Levine_13[["ACCENSE"]] <- (5 * 60) + 7
runtime_Nilsson[["ACCENSE"]] <- (5 * 60) + 39
runtime_Mosmann[["ACCENSE"]] <- (4 * 60) + 17

runtime_Levine_32[["FLOCK"]] <- 30
runtime_Levine_13[["FLOCK"]] <- 7
runtime_Nilsson[["FLOCK"]] <- 6
runtime_Mosmann[["FLOCK"]] <- 65

runtime_Levine_32[["PhenoGraph"]] <- (6 * 60) + 41
runtime_Levine_13[["PhenoGraph"]] <- (4 * 60) + 18
runtime_Nilsson[["PhenoGraph"]] <- (2 * 60) + 1
runtime_Mosmann[["PhenoGraph"]] <- (49 * 60) + 32

runtime_Levine_32[["SWIFT"]] <- 3600 + (7 * 60) + 47
runtime_Levine_13[["SWIFT"]] <- (25 * 60) + 16
runtime_Nilsson[["SWIFT"]] <- (13 * 60) + 5
runtime_Mosmann[["SWIFT"]] <- (41 * 60) + 22



# ==========================================
# Load runtime results for all other methods
# ==========================================

RUNTIME_DIR <- "../results/runtime"


other_names_Levine_32 <- c("DensVM", "flowMeans", "FlowSOM", "FlowSOM_meta", "immunoClust", 
                           "immunoClust_all", "kmeans", "Rclusterpp", "SamSPECTRAL")

other_names_Levine_13 <- c("DensVM", "flowMeans", "FlowSOM", "FlowSOM_meta", "immunoClust", 
                           "immunoClust_all", "kmeans", "Rclusterpp", "SamSPECTRAL")

other_names_Nilsson <- c("DensVM", "flowMeans", "FlowSOM", "FlowSOM_meta", "immunoClust", 
                         "immunoClust_all", "kmeans", "Rclusterpp", "SamSPECTRAL")
# no Rclusterpp
other_names_Mosmann <- c("DensVM", "flowMeans", "FlowSOM", "FlowSOM_meta", "immunoClust", 
                         "immunoClust_all", "kmeans", "SamSPECTRAL")


# load runtime results files

for (i in 1:length(other_names_Levine_32)) {
  file_i <- paste0(RUNTIME_DIR, "/runtime_", other_names_Levine_32[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  runtime_i <- data_i["Levine_2015_marrow_32", "runtime"]
  runtime_Levine_32[[other_names_Levine_32[i]]] <- runtime_i
}

for (i in 1:length(other_names_Levine_13)) {
  file_i <- paste0(RUNTIME_DIR, "/runtime_", other_names_Levine_13[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  runtime_i <- data_i["Levine_2015_marrow_13", "runtime"]
  runtime_Levine_13[[other_names_Levine_13[i]]] <- runtime_i
}

for (i in 1:length(other_names_Nilsson)) {
  file_i <- paste0(RUNTIME_DIR, "/runtime_", other_names_Nilsson[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  runtime_i <- data_i["Nilsson_2013_HSC", "runtime"]
  runtime_Nilsson[[other_names_Nilsson[i]]] <- runtime_i
}

for (i in 1:length(other_names_Mosmann)) {
  file_i <- paste0(RUNTIME_DIR, "/runtime_", other_names_Mosmann[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  runtime_i <- data_i["Mosmann_2014_activ", "runtime"]
  runtime_Mosmann[[other_names_Mosmann[i]]] <- runtime_i
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


