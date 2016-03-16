#########################################################################################
# R script to run DensVM (cytofkit R/Bioconductor package)
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)
library(cytofkit)



#################
### LOAD DATA ###
#################

# working directory (need to reset after running DensVM)

INIT_DIR <- getwd()

# use non-transformed data files, since DensVM will transform automatically

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransform.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13_notransform.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC_notransform.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ_notransform.fcs")

data_Levine_32 <- flowCore::read.FCS(file_Levine_32, transformation = FALSE)
data_Levine_13 <- flowCore::read.FCS(file_Levine_13, transformation = FALSE)
data_Nilsson <- flowCore::read.FCS(file_Nilsson, transformation = FALSE)
data_Mosmann <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)

head(data_Levine_32)
head(data_Levine_13)
head(data_Nilsson)
head(data_Mosmann)

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)


# note DensVM also requires a copy of the input FCS file in the output directory

out_dir_Levine_32 <- "../results_auto/DensVM/Levine_2015_marrow_32"
out_dir_Levine_13 <- "../results_auto/DensVM/Levine_2015_marrow_13"
out_dir_Nilsson <- "../results_auto/DensVM/Nilsson_2013_HSC"
out_dir_Mosmann <- "../results_auto/DensVM/Mosmann_2014_activ"


# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)

# select parameters (protein names) for DensVM

para_Levine_32 <- colnames(data_Levine_32)[marker_cols_Levine_32]
para_Levine_13 <- colnames(data_Levine_13)[marker_cols_Levine_13]
para_Nilsson <- colnames(data_Nilsson)[marker_cols_Nilsson]
para_Mosmann <- colnames(data_Mosmann)[marker_cols_Mosmann]

length(para_Levine_32)
length(para_Levine_13)
length(para_Nilsson)
length(para_Mosmann)

# number of points to subsample

n_sub_Levine_32 <- 20000
n_sub_Levine_13 <- 20000
n_sub_Nilsson <- 20000
n_sub_Mosmann <- 20000




################################################
### Run DensVM: automatic number of clusters ###
################################################

# run DensVM with automatic selection of number of clusters (manual selection is not available)

# use main function "cytof_tsne_densvm"
# graphical version can also be launched with "cytof_tsne_densvm_GUI()"

# results will be saved to files in the output directories

set.seed(123)
runtime_Levine_32 <- system.time({
  cytofkit::cytof_tsne_densvm(fcsFile = file_Levine_32, resDir = out_dir_Levine_32, 
                              lgclMethod = "auto", para = para_Levine_32, 
                              fixedNum = n_sub_Levine_32, verbose = TRUE)
})
setwd(INIT_DIR)  # reset working directory


set.seed(123)
runtime_Levine_13 <- system.time({
  cytofkit::cytof_tsne_densvm(fcsFile = file_Levine_13, resDir = out_dir_Levine_13, 
                              lgclMethod = "auto", para = para_Levine_13, 
                              fixedNum = n_sub_Levine_13, verbose = TRUE)
})
setwd(INIT_DIR)  # reset working directory


set.seed(123)
runtime_Nilsson <- system.time({
  cytofkit::cytof_tsne_densvm(fcsFile = file_Nilsson, resDir = out_dir_Nilsson, 
                              lgclMethod = "auto", para = para_Nilsson, 
                              fixedNum = n_sub_Nilsson, verbose = TRUE)
})
setwd(INIT_DIR)  # reset working directory


set.seed(123)
runtime_Mosmann <- system.time({
  cytofkit::cytof_tsne_densvm(fcsFile = file_Mosmann, resDir = out_dir_Mosmann, 
                              lgclMethod = "auto", para = para_Mosmann, 
                              fixedNum = n_sub_Mosmann, verbose = TRUE)
})
setwd(INIT_DIR)  # reset working directory


# save runtime

runtime_DensVM <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_DensVM, 
            file = "../results_auto/runtime/runtime_DensVM.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_auto/session_info/session_info_DensVM.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_auto/RData_files/results_DensVM.RData")


