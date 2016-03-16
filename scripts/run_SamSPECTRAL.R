#########################################################################################
# R script to run SamSPECTRAL
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)
library(SamSPECTRAL)



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")

data_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_Levine_32, transformation = FALSE))
data_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_Levine_13, transformation = FALSE))
data_Nilsson <- flowCore::exprs(flowCore::read.FCS(file_Nilsson, transformation = FALSE))
data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

head(data_Levine_32)
head(data_Levine_13)
head(data_Nilsson)
head(data_Mosmann)

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)


# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)


# subset data

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
data_Nilsson <- data_Nilsson[, marker_cols_Nilsson]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)




#####################################################
### Run SamSPECTRAL: automatic number of clusters ###
#####################################################

# run SamSPECTRAL with selection of number of clusters via arguments normal.sigma and separation.factor
# (values from Levine et al. 2015)

set.seed(123)
runtime_Levine_32_auto <- system.time({
  out_SamSPECTRAL_Levine_32_auto <- SamSPECTRAL(data_Levine_32, normal.sigma = 100, separation.factor = 1)
})


set.seed(123)
runtime_Levine_13_auto <- system.time({
  out_SamSPECTRAL_Levine_13_auto <- SamSPECTRAL(data_Levine_13, normal.sigma = 100, separation.factor = 1)
})


set.seed(123)
runtime_Nilsson_auto <- system.time({
  out_SamSPECTRAL_Nilsson_auto <- SamSPECTRAL(data_Nilsson, normal.sigma = 100, separation.factor = 1)
})


set.seed(123)
runtime_Mosmann_auto <- system.time({
  out_SamSPECTRAL_Mosmann_auto <- SamSPECTRAL(data_Mosmann, normal.sigma = 100, separation.factor = 1)
})


# extract cluster labels

clus_SamSPECTRAL_Levine_32_auto <- out_SamSPECTRAL_Levine_32_auto
clus_SamSPECTRAL_Levine_13_auto <- out_SamSPECTRAL_Levine_13_auto
clus_SamSPECTRAL_Nilsson_auto <- out_SamSPECTRAL_Nilsson_auto
clus_SamSPECTRAL_Mosmann_auto <- out_SamSPECTRAL_Mosmann_auto

length(clus_SamSPECTRAL_Levine_32_auto)
length(clus_SamSPECTRAL_Levine_13_auto)
length(clus_SamSPECTRAL_Nilsson_auto)
length(clus_SamSPECTRAL_Mosmann_auto)


# cluster sizes and number of clusters

table(clus_SamSPECTRAL_Levine_32_auto)
table(clus_SamSPECTRAL_Levine_13_auto)
table(clus_SamSPECTRAL_Nilsson_auto)
table(clus_SamSPECTRAL_Mosmann_auto)

length(table(clus_SamSPECTRAL_Levine_32_auto))
length(table(clus_SamSPECTRAL_Levine_13_auto))
length(table(clus_SamSPECTRAL_Nilsson_auto))
length(table(clus_SamSPECTRAL_Mosmann_auto))


# save cluster labels

res_SamSPECTRAL_Levine_32_auto <- data.frame(label = clus_SamSPECTRAL_Levine_32_auto)
res_SamSPECTRAL_Levine_13_auto <- data.frame(label = clus_SamSPECTRAL_Levine_13_auto)
res_SamSPECTRAL_Nilsson_auto <- data.frame(label = clus_SamSPECTRAL_Nilsson_auto)
res_SamSPECTRAL_Mosmann_auto <- data.frame(label = clus_SamSPECTRAL_Mosmann_auto)


write.table(res_SamSPECTRAL_Levine_32_auto, 
            file = "../results_auto/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Levine_13_auto, 
            file = "../results_auto/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Nilsson_auto, 
            file = "../results_auto/SamSPECTRAL/SamSPECTRAL_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Mosmann_auto, 
            file = "../results_auto/SamSPECTRAL/SamSPECTRAL_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_SamSPECTRAL_auto <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32_auto["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13_auto["elapsed"], 
  Nilsson_2013_HSC = runtime_Nilsson_auto["elapsed"], 
  Mosmann_2014_activ = runtime_Mosmann_auto["elapsed"], 
  row.names = "runtime"))

write.table(runtime_SamSPECTRAL_auto, 
            file = "../results_auto/runtime/runtime_SamSPECTRAL.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_auto/session_info/session_info_SamSPECTRAL.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_auto/RData_files/results_SamSPECTRAL.RData")




#############################################################
### Run SamSPECTRAL: manually selected number of clusters ###
#############################################################

# run SamSPECTRAL with manual selection of number of clusters


# number of clusters
k_Levine_32 <- 40
k_Levine_13 <- 40
k_Nilsson <- 40
k_Mosmann <- 40


set.seed(123)
runtime_Levine_32_manual <- system.time({
  out_SamSPECTRAL_Levine_32_manual <- SamSPECTRAL(data_Levine_32, normal.sigma = 100, separation.factor = 1, 
                                                  number.of.clusters = k_Levine_32)
})


set.seed(123)
runtime_Levine_13_manual <- system.time({
  out_SamSPECTRAL_Levine_13_manual <- SamSPECTRAL(data_Levine_13, normal.sigma = 100, separation.factor = 1, 
                                                  number.of.clusters = k_Levine_13)
})


set.seed(123)
runtime_Nilsson_manual <- system.time({
  out_SamSPECTRAL_Nilsson_manual <- SamSPECTRAL(data_Nilsson, normal.sigma = 100, separation.factor = 1, 
                                                number.of.clusters = k_Nilsson)
})


set.seed(123)
runtime_Mosmann_manual <- system.time({
  out_SamSPECTRAL_Mosmann_manual <- SamSPECTRAL(data_Mosmann, normal.sigma = 100, separation.factor = 1, 
                                                number.of.clusters = k_Mosmann)
})


# extract cluster labels

clus_SamSPECTRAL_Levine_32_manual <- out_SamSPECTRAL_Levine_32_manual
clus_SamSPECTRAL_Levine_13_manual <- out_SamSPECTRAL_Levine_13_manual
clus_SamSPECTRAL_Nilsson_manual <- out_SamSPECTRAL_Nilsson_manual
clus_SamSPECTRAL_Mosmann_manual <- out_SamSPECTRAL_Mosmann_manual

length(clus_SamSPECTRAL_Levine_32_manual)
length(clus_SamSPECTRAL_Levine_13_manual)
length(clus_SamSPECTRAL_Nilsson_manual)
length(clus_SamSPECTRAL_Mosmann_manual)


# cluster sizes and number of clusters

table(clus_SamSPECTRAL_Levine_32_manual)
table(clus_SamSPECTRAL_Levine_13_manual)
table(clus_SamSPECTRAL_Nilsson_manual)
table(clus_SamSPECTRAL_Mosmann_manual)

length(table(clus_SamSPECTRAL_Levine_32_manual))
length(table(clus_SamSPECTRAL_Levine_13_manual))
length(table(clus_SamSPECTRAL_Nilsson_manual))
length(table(clus_SamSPECTRAL_Mosmann_manual))


# save cluster labels

res_SamSPECTRAL_Levine_32_manual <- data.frame(label = clus_SamSPECTRAL_Levine_32_manual)
res_SamSPECTRAL_Levine_13_manual <- data.frame(label = clus_SamSPECTRAL_Levine_13_manual)
res_SamSPECTRAL_Nilsson_manual <- data.frame(label = clus_SamSPECTRAL_Nilsson_manual)
res_SamSPECTRAL_Mosmann_manual <- data.frame(label = clus_SamSPECTRAL_Mosmann_manual)


write.table(res_SamSPECTRAL_Levine_32_manual, 
            file = "../results_manual/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Levine_13_manual, 
            file = "../results_manual/SamSPECTRAL/SamSPECTRAL_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Nilsson_manual, 
            file = "../results_manual/SamSPECTRAL/SamSPECTRAL_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_SamSPECTRAL_Mosmann_manual, 
            file = "../results_manual/SamSPECTRAL/SamSPECTRAL_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_SamSPECTRAL_manual <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32_manual["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13_manual["elapsed"], 
  Nilsson_2013_HSC = runtime_Nilsson_manual["elapsed"], 
  Mosmann_2014_activ = runtime_Mosmann_manual["elapsed"], 
  row.names = "runtime"))

write.table(runtime_SamSPECTRAL_manual, 
            file = "../results_manual/runtime/runtime_SamSPECTRAL.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_manual/session_info/session_info_SamSPECTRAL.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_manual/RData_files/results_SamSPECTRAL.RData")


