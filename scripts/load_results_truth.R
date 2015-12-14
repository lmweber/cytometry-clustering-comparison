#########################################################################################
# R script to load results for manually gated populations (truth)
#
# Lukas M. Weber, December 2015
#########################################################################################


library(flowCore)



#################################################
### load true (manually gated) cluster labels ###
#################################################

DATA_DIR <- "../../benchmark_data_sets/"

file_truth_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_truth_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_truth_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_truth_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")

data_truth_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_32, transformation = FALSE))
data_truth_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_truth_Levine_13, transformation = FALSE))
data_truth_Nilsson <- flowCore::exprs(flowCore::read.FCS(file_truth_Nilsson, transformation = FALSE))
data_truth_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_truth_Mosmann, transformation = FALSE))

dim(data_truth_Levine_32)
dim(data_truth_Levine_13)
dim(data_truth_Nilsson)
dim(data_truth_Mosmann)


# extract cluster labels

clus_truth_Levine_32 <- data_truth_Levine_32[, "label"]
clus_truth_Levine_13 <- data_truth_Levine_13[, "label"]
clus_truth_Nilsson <- data_truth_Nilsson[, "label"]
clus_truth_Mosmann <- data_truth_Mosmann[, "label"]

length(clus_truth_Levine_32)
length(clus_truth_Levine_13)
length(clus_truth_Nilsson)
length(clus_truth_Mosmann)


# cluster sizes and number of clusters

tbl_truth_Levine_32 <- table(clus_truth_Levine_32)
tbl_truth_Levine_13 <- table(clus_truth_Levine_13)
tbl_truth_Nilsson <- table(clus_truth_Nilsson)
tbl_truth_Mosmann <- table(clus_truth_Mosmann)

tbl_truth_Levine_32
tbl_truth_Levine_13
tbl_truth_Nilsson
tbl_truth_Mosmann

length(tbl_truth_Levine_32)  # 14 populations
length(tbl_truth_Levine_13)  # 24 populations
length(tbl_truth_Nilsson)  # 1 = rare population of interest, 0 = all others
length(tbl_truth_Mosmann)  # 1 = rare population of interest, 0 = all others


