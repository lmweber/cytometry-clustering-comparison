#########################################################################################
# R script to run SamSPECTRAL
#
# Lukas M. Weber, December 2015
#########################################################################################


library(flowCore)
library(SamSPECTRAL)



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare.fcs")

data_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_Levine_32, transformation = FALSE))
data_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_Levine_13, transformation = FALSE))
data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

head(data_Levine_32)
head(data_Levine_13)
head(data_Mosmann)

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Mosmann)

# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Mosmann)

# subset data

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Mosmann)



#######################
### Run SamSPECTRAL ###
#######################

# run SamSPECTRAL

set.seed(123)
runtime_Levine_32 <- system.time({
  out_SamSPECTRAL_Levine_32 <- SamSPECTRAL(data_Levine_32, normal.sigma = 100, separation.factor = 1)
})

set.seed(123)
runtime_Levine_13 <- system.time({
  out_SamSPECTRAL_Levine_13 <- SamSPECTRAL(data_Levine_13, normal.sigma = 100, separation.factor = 1)
})

set.seed(123)
runtime_Mosmann <- system.time({
  out_SamSPECTRAL_Mosmann <- SamSPECTRAL(data_Mosmann, normal.sigma = 100, separation.factor = 1)
})


# extract cluster labels

clus_SamSPECTRAL_Levine_32 <- out_SamSPECTRAL_Levine_32
clus_SamSPECTRAL_Levine_13 <- out_SamSPECTRAL_Levine_13
clus_SamSPECTRAL_Mosmann <- out_SamSPECTRAL_Mosmann

length(clus_SamSPECTRAL_Levine_32)
length(clus_SamSPECTRAL_Levine_13)
length(clus_SamSPECTRAL_Mosmann)


# cluster sizes and number of clusters

table(clus_SamSPECTRAL_Levine_32)
table(clus_SamSPECTRAL_Levine_13)
table(clus_SamSPECTRAL_Mosmann)

length(table(clus_SamSPECTRAL_Levine_32))
length(table(clus_SamSPECTRAL_Levine_13))
length(table(clus_SamSPECTRAL_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_SamSPECTRAL_Levine_32 <- data.frame(label = clus_SamSPECTRAL_Levine_32)
res_SamSPECTRAL_Levine_13 <- data.frame(label = clus_SamSPECTRAL_Levine_13)
res_SamSPECTRAL_Mosmann <- data.frame(label = clus_SamSPECTRAL_Mosmann)

write.table(res_SamSPECTRAL_Levine_32, 
            file = "../results/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Levine_13, 
            file = "../results/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Mosmann, 
            file = "../results/SamSPECTRAL/SamSPECTRAL_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_SamSPECTRAL <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13["elapsed"], 
  Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_SamSPECTRAL, file = "../results/runtime/runtime_SamSPECTRAL.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/SamSPECTRAL_session_info.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/SamSPECTRAL_results.RData")


