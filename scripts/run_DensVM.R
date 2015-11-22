#########################################################################################
# R script to run DensVM (cytofkit R/Bioconductor package)
#
# Lukas M. Weber, November 2015
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

file_Levine <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransf.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare_notransf.fcs")

data_Levine <- flowCore::read.FCS(file_Levine, transformation = FALSE)
data_Mosmann <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)

head(data_Levine)
head(data_Mosmann)

dim(data_Levine)
dim(data_Mosmann)

# note DensVM also requires a copy of the input FCS file in the output directory

out_dir_Levine <- "../results/DensVM/Levine_2015_marrow_32"
out_dir_Mosmann <- "../results/DensVM/Mosmann_2014_rare"

# indices of protein marker columns

marker_cols_Levine <- 5:36
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine)
length(marker_cols_Mosmann)

# select parameters (protein names) for DensVM

para_Levine <- colnames(data_Levine)[marker_cols_Levine]
para_Mosmann <- colnames(data_Mosmann)[marker_cols_Mosmann]

length(para_Levine)
length(para_Mosmann)

# number of points to subsample

n_sub_Levine <- 20000
n_sub_Mosmann <- 20000



##################
### Run DensVM ###
##################

# use main function "cytof_tsne_densvm"
# graphical version can also be launched with "cytof_tsne_densvm_GUI()"

# results will be saved to files in the output directories

set.seed(123)
runtime_Levine <- system.time({
  cytofkit::cytof_tsne_densvm(fcsFile = file_Levine, resDir = out_dir_Levine, 
                              para = para_Levine, fixedNum = n_sub_Levine, verbose = TRUE)
})
setwd(INIT_DIR)  # reset working directory


set.seed(123)
runtime_Mosmann <- system.time({
  cytofkit::cytof_tsne_densvm(fcsFile = file_Mosmann, resDir = out_dir_Mosmann, 
                              para = para_Mosmann, fixedNum = n_sub_Mosmann, verbose = TRUE)
})
setwd(INIT_DIR)  # reset working directory



####################
### SAVE RESULTS ###
####################

# save runtime

runtime_DensVM <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine["elapsed"], 
  Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_DensVM, file = "../results/runtime/runtime_DensVM.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/DensVM_session_info.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/DensVM_results.RData")


