#########################################################################################
# R script to run FlowSOM and FlowSOM_meta
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)
library(FlowSOM)



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ.fcs")


# FlowSOM requires input data as flowFrame objects

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


# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- c(7:9, 11:21)

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)




#################################################
### Run FlowSOM: automatic number of clusters ###
#################################################

# run FlowSOM with default number of clusters for all data sets (10x10 grid or 100 clusters)

set.seed(123)
runtime_FlowSOM_Levine_32 <- system.time({
  fSOM_Levine_32 <- FlowSOM::ReadInput(data_Levine_32, transform = FALSE, scale = FALSE)
  fSOM_Levine_32 <- FlowSOM::BuildSOM(fSOM_Levine_32, colsToUse = marker_cols_Levine_32)
  fSOM_Levine_32 <- FlowSOM::BuildMST(fSOM_Levine_32)
})


set.seed(123)
runtime_FlowSOM_Levine_13 <- system.time({
  fSOM_Levine_13 <- FlowSOM::ReadInput(data_Levine_13, transform = FALSE, scale = FALSE)
  fSOM_Levine_13 <- FlowSOM::BuildSOM(fSOM_Levine_13, colsToUse = marker_cols_Levine_13)
  fSOM_Levine_13 <- FlowSOM::BuildMST(fSOM_Levine_13)
})


set.seed(100)
runtime_FlowSOM_Nilsson <- system.time({
  fSOM_Nilsson <- FlowSOM::ReadInput(data_Nilsson, transform = FALSE, scale = FALSE)
  fSOM_Nilsson <- FlowSOM::BuildSOM(fSOM_Nilsson, colsToUse = marker_cols_Nilsson)
  fSOM_Nilsson <- FlowSOM::BuildMST(fSOM_Nilsson)
})


set.seed(123)
runtime_FlowSOM_Mosmann <- system.time({
  fSOM_Mosmann <- FlowSOM::ReadInput(data_Mosmann, transform = FALSE, scale = FALSE)
  fSOM_Mosmann <- FlowSOM::BuildSOM(fSOM_Mosmann, colsToUse = marker_cols_Mosmann)
  fSOM_Mosmann <- FlowSOM::BuildMST(fSOM_Mosmann)
})


# plots

FlowSOM::PlotStars(fSOM_Levine_32)
FlowSOM::PlotStars(fSOM_Levine_13)
FlowSOM::PlotStars(fSOM_Nilsson)
FlowSOM::PlotStars(fSOM_Mosmann)


# extract cluster labels

str(fSOM_Levine_32$map)

head(fSOM_Levine_32$map$mapping)
dim(fSOM_Levine_32$map$mapping)


clus_FlowSOM_Levine_32 <- fSOM_Levine_32$map$mapping[, 1]
clus_FlowSOM_Levine_13 <- fSOM_Levine_13$map$mapping[, 1]
clus_FlowSOM_Nilsson <- fSOM_Nilsson$map$mapping[, 1]
clus_FlowSOM_Mosmann <- fSOM_Mosmann$map$mapping[, 1]

length(clus_FlowSOM_Levine_32)
length(clus_FlowSOM_Levine_13)
length(clus_FlowSOM_Nilsson)
length(clus_FlowSOM_Mosmann)


# cluster sizes and number of clusters

table(clus_FlowSOM_Levine_32)
table(clus_FlowSOM_Levine_13)
table(clus_FlowSOM_Nilsson)
table(clus_FlowSOM_Mosmann)

length(table(clus_FlowSOM_Levine_32))
length(table(clus_FlowSOM_Levine_13))
length(table(clus_FlowSOM_Nilsson))
length(table(clus_FlowSOM_Mosmann))


# save cluster labels

res_FlowSOM_Levine_32 <- data.frame(label = clus_FlowSOM_Levine_32)
res_FlowSOM_Levine_13 <- data.frame(label = clus_FlowSOM_Levine_13)
res_FlowSOM_Nilsson <- data.frame(label = clus_FlowSOM_Nilsson)
res_FlowSOM_Mosmann <- data.frame(label = clus_FlowSOM_Mosmann)

write.table(res_FlowSOM_Levine_32, 
            file = "../results_auto/FlowSOM/FlowSOM_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Levine_13, 
            file = "../results_auto/FlowSOM/FlowSOM_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Nilsson, 
            file = "../results_auto/FlowSOM/FlowSOM_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Mosmann, 
            file = "../results_auto/FlowSOM/FlowSOM_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_FlowSOM <- t(data.frame(
  Levine_2015_marrow_32 = runtime_FlowSOM_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_FlowSOM_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_FlowSOM_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_FlowSOM_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_FlowSOM, file = "../results_auto/runtime/runtime_FlowSOM.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_auto/session_info/session_info_FlowSOM.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_auto/RData_files/results_FlowSOM.RData")




######################################################
### Run FlowSOM_meta: automatic number of clusters ###
######################################################

# run metaclustering step (i.e. FlowSOM_meta) with automatic selection of number of clusters

set.seed(123)
runtime_FlowSOM_meta_Levine_32_auto <- system.time({
  meta_clustering_Levine_32_auto <- FlowSOM::MetaClustering(fSOM_Levine_32$map$codes, method = "metaClustering_consensus")
})

set.seed(123)
runtime_FlowSOM_meta_Levine_13_auto <- system.time({
  meta_clustering_Levine_13_auto <- FlowSOM::MetaClustering(fSOM_Levine_13$map$codes, method = "metaClustering_consensus")
})

set.seed(100)
runtime_FlowSOM_meta_Nilsson_auto <- system.time({
  meta_clustering_Nilsson_auto <- FlowSOM::MetaClustering(fSOM_Nilsson$map$codes, method = "metaClustering_consensus")
})

set.seed(123)
runtime_FlowSOM_meta_Mosmann_auto <- system.time({
  meta_clustering_Mosmann_auto <- FlowSOM::MetaClustering(fSOM_Mosmann$map$codes, method = "metaClustering_consensus")
})


# combine runtime

runtime_FlowSOM_meta_Levine_32_auto <- runtime_FlowSOM_Levine_32 + runtime_FlowSOM_meta_Levine_32_auto
runtime_FlowSOM_meta_Levine_13_auto <- runtime_FlowSOM_Levine_13 + runtime_FlowSOM_meta_Levine_13_auto
runtime_FlowSOM_meta_Nilsson_auto <- runtime_FlowSOM_Nilsson + runtime_FlowSOM_meta_Nilsson_auto
runtime_FlowSOM_meta_Mosmann_auto <- runtime_FlowSOM_Mosmann + runtime_FlowSOM_meta_Mosmann_auto


# extract cluster labels

meta_clustering_Levine_32_auto
meta_clustering_Levine_13_auto
meta_clustering_Nilsson_auto
meta_clustering_Mosmann_auto

clus_FlowSOM_meta_Levine_32_auto <- meta_clustering_Levine_32_auto[fSOM_Levine_32$map$mapping[, 1]]
clus_FlowSOM_meta_Levine_13_auto <- meta_clustering_Levine_13_auto[fSOM_Levine_13$map$mapping[, 1]]
clus_FlowSOM_meta_Nilsson_auto <- meta_clustering_Nilsson_auto[fSOM_Nilsson$map$mapping[, 1]]
clus_FlowSOM_meta_Mosmann_auto <- meta_clustering_Mosmann_auto[fSOM_Mosmann$map$mapping[, 1]]

length(clus_FlowSOM_meta_Levine_32_auto)
length(clus_FlowSOM_meta_Levine_13_auto)
length(clus_FlowSOM_meta_Nilsson_auto)
length(clus_FlowSOM_meta_Mosmann_auto)


# cluster sizes and number of clusters

table(clus_FlowSOM_meta_Levine_32_auto)
table(clus_FlowSOM_meta_Levine_13_auto)
table(clus_FlowSOM_meta_Nilsson_auto)
table(clus_FlowSOM_meta_Mosmann_auto)

length(table(clus_FlowSOM_meta_Levine_32_auto))
length(table(clus_FlowSOM_meta_Levine_13_auto))
length(table(clus_FlowSOM_meta_Nilsson_auto))
length(table(clus_FlowSOM_meta_Mosmann_auto))


# save cluster labels

res_FlowSOM_meta_Levine_32_auto <- data.frame(label = clus_FlowSOM_meta_Levine_32_auto)
res_FlowSOM_meta_Levine_13_auto <- data.frame(label = clus_FlowSOM_meta_Levine_13_auto)
res_FlowSOM_meta_Nilsson_auto <- data.frame(label = clus_FlowSOM_meta_Nilsson_auto)
res_FlowSOM_meta_Mosmann_auto <- data.frame(label = clus_FlowSOM_meta_Mosmann_auto)

write.table(res_FlowSOM_meta_Levine_32_auto, 
            file = "../results_auto/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Levine_13_auto, 
            file = "../results_auto/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Nilsson_auto, 
            file = "../results_auto/FlowSOM_meta/FlowSOM_meta_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Mosmann_auto, 
            file = "../results_auto/FlowSOM_meta/FlowSOM_meta_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_FlowSOM_meta_auto <- t(data.frame(
  Levine_2015_marrow_32 = runtime_FlowSOM_meta_Levine_32_auto["elapsed"], 
  Levine_2015_marrow_13 = runtime_FlowSOM_meta_Levine_13_auto["elapsed"], 
  Nilsson_2013_HSC = runtime_FlowSOM_meta_Nilsson_auto["elapsed"], 
  Mosmann_2014_activ = runtime_FlowSOM_meta_Mosmann_auto["elapsed"], 
  row.names = "runtime"))

write.table(runtime_FlowSOM_meta_auto, 
            file = "../results_auto/runtime/runtime_FlowSOM_meta.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_auto/session_info/session_info_FlowSOM_meta.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_auto/RData_files/results_FlowSOM_meta.RData")




#########################################################
### Run FlowSOM: manually selected number of clusters ###
#########################################################

# run FlowSOM with manually selected number of clusters (e.g. 10x10 or 20x20 grid)


# grid size (e.g. 10x10 or 20x20 grid, i.e. 100 or 400 clusters)
grid_Levine_32 <- 10
grid_Levine_13 <- 10
grid_Nilsson <- 10
grid_Mosmann <- 20


set.seed(123)
runtime_FlowSOM_Levine_32 <- system.time({
  fSOM_Levine_32 <- FlowSOM::ReadInput(data_Levine_32, transform = FALSE, scale = FALSE)
  fSOM_Levine_32 <- FlowSOM::BuildSOM(fSOM_Levine_32, colsToUse = marker_cols_Levine_32, 
                                      xdim = grid_Levine_32, ydim = grid_Levine_32)
  fSOM_Levine_32 <- FlowSOM::BuildMST(fSOM_Levine_32)
})


set.seed(123)
runtime_FlowSOM_Levine_13 <- system.time({
  fSOM_Levine_13 <- FlowSOM::ReadInput(data_Levine_13, transform = FALSE, scale = FALSE)
  fSOM_Levine_13 <- FlowSOM::BuildSOM(fSOM_Levine_13, colsToUse = marker_cols_Levine_13, 
                                      xdim = grid_Levine_13, ydim = grid_Levine_13)
  fSOM_Levine_13 <- FlowSOM::BuildMST(fSOM_Levine_13)
})


set.seed(100)
runtime_FlowSOM_Nilsson <- system.time({
  fSOM_Nilsson <- FlowSOM::ReadInput(data_Nilsson, transform = FALSE, scale = FALSE)
  fSOM_Nilsson <- FlowSOM::BuildSOM(fSOM_Nilsson, colsToUse = marker_cols_Nilsson, 
                                    xdim = grid_Nilsson, ydim = grid_Nilsson)
  fSOM_Nilsson <- FlowSOM::BuildMST(fSOM_Nilsson)
})


set.seed(123)
runtime_FlowSOM_Mosmann <- system.time({
  fSOM_Mosmann <- FlowSOM::ReadInput(data_Mosmann, transform = FALSE, scale = FALSE)
  fSOM_Mosmann <- FlowSOM::BuildSOM(fSOM_Mosmann, colsToUse = marker_cols_Mosmann, 
                                    xdim = grid_Mosmann, ydim = grid_Mosmann)
  fSOM_Mosmann <- FlowSOM::BuildMST(fSOM_Mosmann)
})


# plots

FlowSOM::PlotStars(fSOM_Levine_32)
FlowSOM::PlotStars(fSOM_Levine_13)
FlowSOM::PlotStars(fSOM_Nilsson)
FlowSOM::PlotStars(fSOM_Mosmann)


# extract cluster labels

str(fSOM_Levine_32$map)

head(fSOM_Levine_32$map$mapping)
dim(fSOM_Levine_32$map$mapping)


clus_FlowSOM_Levine_32 <- fSOM_Levine_32$map$mapping[, 1]
clus_FlowSOM_Levine_13 <- fSOM_Levine_13$map$mapping[, 1]
clus_FlowSOM_Nilsson <- fSOM_Nilsson$map$mapping[, 1]
clus_FlowSOM_Mosmann <- fSOM_Mosmann$map$mapping[, 1]

length(clus_FlowSOM_Levine_32)
length(clus_FlowSOM_Levine_13)
length(clus_FlowSOM_Nilsson)
length(clus_FlowSOM_Mosmann)


# cluster sizes and number of clusters

table(clus_FlowSOM_Levine_32)
table(clus_FlowSOM_Levine_13)
table(clus_FlowSOM_Nilsson)
table(clus_FlowSOM_Mosmann)

length(table(clus_FlowSOM_Levine_32))
length(table(clus_FlowSOM_Levine_13))
length(table(clus_FlowSOM_Nilsson))
length(table(clus_FlowSOM_Mosmann))


# save cluster labels

res_FlowSOM_Levine_32 <- data.frame(label = clus_FlowSOM_Levine_32)
res_FlowSOM_Levine_13 <- data.frame(label = clus_FlowSOM_Levine_13)
res_FlowSOM_Nilsson <- data.frame(label = clus_FlowSOM_Nilsson)
res_FlowSOM_Mosmann <- data.frame(label = clus_FlowSOM_Mosmann)

write.table(res_FlowSOM_Levine_32, 
            file = "../results_manual/FlowSOM/FlowSOM_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Levine_13, 
            file = "../results_manual/FlowSOM/FlowSOM_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Nilsson, 
            file = "../results_manual/FlowSOM/FlowSOM_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Mosmann, 
            file = "../results_manual/FlowSOM/FlowSOM_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_FlowSOM <- t(data.frame(
  Levine_2015_marrow_32 = runtime_FlowSOM_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_FlowSOM_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_FlowSOM_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_FlowSOM_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_FlowSOM, file = "../results_manual/runtime/runtime_FlowSOM.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_manual/session_info/session_info_FlowSOM.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_manual/RData_files/results_FlowSOM.RData")




##############################################################
### Run FlowSOM_meta: manually selected number of clusters ###
##############################################################

# run metaclustering step (i.e. FlowSOM_meta) with manual selection of number of clusters


# number of clusters
k_Levine_32 <- 40
k_Levine_13 <- 40
k_Nilsson <- 40
k_Mosmann <- 40


set.seed(123)
runtime_FlowSOM_meta_Levine_32_manual <- system.time({
  meta_clustering_Levine_32_manual <- FlowSOM::metaClustering_consensus(fSOM_Levine_32$map$codes, k = k_Levine_32)
})

set.seed(123)
runtime_FlowSOM_meta_Levine_13_manual <- system.time({
  meta_clustering_Levine_13_manual <- FlowSOM::metaClustering_consensus(fSOM_Levine_13$map$codes, k = k_Levine_13)
})

set.seed(100)
runtime_FlowSOM_meta_Nilsson_manual <- system.time({
  meta_clustering_Nilsson_manual <- FlowSOM::metaClustering_consensus(fSOM_Nilsson$map$codes, k = k_Nilsson)
})

set.seed(123)
runtime_FlowSOM_meta_Mosmann_manual <- system.time({
  meta_clustering_Mosmann_manual <- FlowSOM::metaClustering_consensus(fSOM_Mosmann$map$codes, k = k_Mosmann)
})


# combine runtime

runtime_FlowSOM_meta_Levine_32_manual <- runtime_FlowSOM_Levine_32 + runtime_FlowSOM_meta_Levine_32_manual
runtime_FlowSOM_meta_Levine_13_manual <- runtime_FlowSOM_Levine_13 + runtime_FlowSOM_meta_Levine_13_manual
runtime_FlowSOM_meta_Nilsson_manual <- runtime_FlowSOM_Nilsson + runtime_FlowSOM_meta_Nilsson_manual
runtime_FlowSOM_meta_Mosmann_manual <- runtime_FlowSOM_Mosmann + runtime_FlowSOM_meta_Mosmann_manual


# extract cluster labels

meta_clustering_Levine_32_manual
meta_clustering_Levine_13_manual
meta_clustering_Nilsson_manual
meta_clustering_Mosmann_manual

clus_FlowSOM_meta_Levine_32_manual <- meta_clustering_Levine_32_manual[fSOM_Levine_32$map$mapping[, 1]]
clus_FlowSOM_meta_Levine_13_manual <- meta_clustering_Levine_13_manual[fSOM_Levine_13$map$mapping[, 1]]
clus_FlowSOM_meta_Nilsson_manual <- meta_clustering_Nilsson_manual[fSOM_Nilsson$map$mapping[, 1]]
clus_FlowSOM_meta_Mosmann_manual <- meta_clustering_Mosmann_manual[fSOM_Mosmann$map$mapping[, 1]]

length(clus_FlowSOM_meta_Levine_32_manual)
length(clus_FlowSOM_meta_Levine_13_manual)
length(clus_FlowSOM_meta_Nilsson_manual)
length(clus_FlowSOM_meta_Mosmann_manual)


# cluster sizes and number of clusters

table(clus_FlowSOM_meta_Levine_32_manual)
table(clus_FlowSOM_meta_Levine_13_manual)
table(clus_FlowSOM_meta_Nilsson_manual)
table(clus_FlowSOM_meta_Mosmann_manual)

length(table(clus_FlowSOM_meta_Levine_32_manual))
length(table(clus_FlowSOM_meta_Levine_13_manual))
length(table(clus_FlowSOM_meta_Nilsson_manual))
length(table(clus_FlowSOM_meta_Mosmann_manual))


# save cluster labels

res_FlowSOM_meta_Levine_32_manual <- data.frame(label = clus_FlowSOM_meta_Levine_32_manual)
res_FlowSOM_meta_Levine_13_manual <- data.frame(label = clus_FlowSOM_meta_Levine_13_manual)
res_FlowSOM_meta_Nilsson_manual <- data.frame(label = clus_FlowSOM_meta_Nilsson_manual)
res_FlowSOM_meta_Mosmann_manual <- data.frame(label = clus_FlowSOM_meta_Mosmann_manual)

write.table(res_FlowSOM_meta_Levine_32_manual, 
            file = "../results_manual/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Levine_13_manual, 
            file = "../results_manual/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Nilsson_manual, 
            file = "../results_manual/FlowSOM_meta/FlowSOM_meta_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Mosmann_manual, 
            file = "../results_manual/FlowSOM_meta/FlowSOM_meta_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_FlowSOM_meta_manual <- t(data.frame(
  Levine_2015_marrow_32 = runtime_FlowSOM_meta_Levine_32_manual["elapsed"], 
  Levine_2015_marrow_13 = runtime_FlowSOM_meta_Levine_13_manual["elapsed"], 
  Nilsson_2013_HSC = runtime_FlowSOM_meta_Nilsson_manual["elapsed"], 
  Mosmann_2014_activ = runtime_FlowSOM_meta_Mosmann_manual["elapsed"], 
  row.names = "runtime"))

write.table(runtime_FlowSOM_meta_manual, 
            file = "../results_manual/runtime/runtime_FlowSOM_meta.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_manual/session_info/session_info_FlowSOM_meta.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_manual/RData_files/results_FlowSOM_meta.RData")


