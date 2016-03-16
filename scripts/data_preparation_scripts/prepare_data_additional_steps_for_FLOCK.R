#########################################################################################
# Additional R script to prepare data files for FLOCK
#
# FLOCK requires TXT data files containing protein expression columns only. This script 
# selects the required columns and exports as TXT files.
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)


### run from FLOCK program directory

DATA_DIR <- "../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")

data_Levine_32 <- exprs(read.FCS(file_Levine_32, transformation = FALSE))
data_Levine_13 <- exprs(read.FCS(file_Levine_13, transformation = FALSE))
data_Nilsson <- exprs(read.FCS(file_Nilsson, transformation = FALSE))
data_Mosmann <- exprs(read.FCS(file_Mosmann, transformation = FALSE))

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- c(7:9, 11:21)

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
data_Nilsson <- data_Nilsson[, marker_cols_Nilsson]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

head(data_Levine_32)
head(data_Levine_13)
head(data_Nilsson)
head(data_Mosmann)

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)

write.table(data_Levine_32, file = "Levine_2015_marrow_32_markers_only.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_Levine_13, file = "Levine_2015_marrow_13_markers_only.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_Nilsson, file = "Nilsson_2013_HSC_markers_only.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_Mosmann, file = "Mosmann_2014_activ_markers_only.txt", quote = FALSE, sep = "\t", row.names = FALSE)


