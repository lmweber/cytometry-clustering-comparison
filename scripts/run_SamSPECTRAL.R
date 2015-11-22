#########################################################################################
# R script to run SamSPECTRAL
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)
library(SamSPECTRAL)



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare.fcs")

data_Levine <- flowCore::exprs(flowCore::read.FCS(file_Levine, transformation = FALSE))
data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

head(data_Levine)
head(data_Mosmann)

dim(data_Levine)
dim(data_Mosmann)

# indices of protein marker columns

marker_cols_Levine <- 5:36
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine)
length(marker_cols_Mosmann)

# subset data

data_Levine <- data_Levine[, marker_cols_Levine]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine)
dim(data_Mosmann)



#######################
### Run SamSPECTRAL ###
#######################

# run SamSPECTRAL

set.seed(123)
runtime_Levine <- system.time({
  out_SamSPECTRAL_Levine <- SamSPECTRAL(data_Levine, normal.sigma = 100, separation.factor = 1)
})

set.seed(123)
runtime_Mosmann <- system.time({
  out_SamSPECTRAL_Mosmann <- SamSPECTRAL(data_Mosmann, normal.sigma = 100, separation.factor = 1)
})


# extract cluster labels

clus_SamSPECTRAL_Levine <- out_SamSPECTRAL_Levine
clus_SamSPECTRAL_Mosmann <- out_SamSPECTRAL_Mosmann

length(clus_SamSPECTRAL_Levine)
length(clus_SamSPECTRAL_Mosmann)


# cluster sizes and number of clusters

table(clus_SamSPECTRAL_Levine)
table(clus_SamSPECTRAL_Mosmann)

length(table(clus_SamSPECTRAL_Levine))
length(table(clus_SamSPECTRAL_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_SamSPECTRAL_Levine <- data.frame(label = clus_SamSPECTRAL_Levine)
res_SamSPECTRAL_Mosmann <- data.frame(label = clus_SamSPECTRAL_Mosmann)

write.table(res_SamSPECTRAL_Levine, 
            file = "../results/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Mosmann, 
            file = "../results/SamSPECTRAL/SamSPECTRAL_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_SamSPECTRAL <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine["elapsed"], 
  Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_SamSPECTRAL, file = "../results/runtime/runtime_SamSPECTRAL.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/SamSPECTRAL_session_info.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/SamSPECTRAL_results.RData")


