#########################################################################################
# R script to load runtimes (including densityCut)
#
# Lukas Weber, February 2017
#########################################################################################


DIR_PREVIOUS <- "../../../evaluate_results"
DIR_DENSITYCUT <- getwd()



# ======================================================
# Load previous results for all other clustering methods
# ======================================================

setwd(DIR_PREVIOUS)

source("evaluate_runtime.R")



# ============================
# Load runtimes for densityCut
# ============================

setwd(DIR_DENSITYCUT)


# run script "evaluate_all_methods.R" to load results directories (automatic or manual 
# number of clusters)

#source("evaluate_all_methods.R")  ## takes 20 min


# load runtime results files

file_i <- paste0(RES_DIR_DENSITYCUT, "/../runtimes/runtime_densityCut.txt")
data_i <- read.table(file_i, header = TRUE, sep = "\t")

res_runtime[["Levine_32dim"]][["densityCut"]] <- data_i["Levine_32dim", "runtime"]
res_runtime[["Levine_13dim"]][["densityCut"]] <- data_i["Levine_13dim", "runtime"]
res_runtime[["Samusik_01"]][["densityCut"]]   <- data_i["Samusik_01", "runtime"]
res_runtime[["Samusik_all"]][["densityCut"]]  <- data_i["Samusik_all", "runtime"]
res_runtime[["Nilsson_rare"]][["densityCut"]] <- data_i["Nilsson_rare", "runtime"]
res_runtime[["Mosmann_rare"]][["densityCut"]] <- data_i["Mosmann_rare", "runtime"]
res_runtime[["FlowCAP_ND"]][["densityCut"]]   <- data_i["FlowCAP_ND", "runtime"]
res_runtime[["FlowCAP_WNV"]][["densityCut"]]  <- data_i["FlowCAP_WNV", "runtime"]



# ===========
# Subsampling
# ===========

# skip this step; no subsampling required for densityCut (see report)

# see script for other clustering methods ("evaluate_runtime.R")



# ========================
# Multiple processor cores
# ========================

# skip this step; multiple processor cores not required for densityCut (see report)

# see script for other clustering methods ("evaluate_runtime.R")


