#########################################################################################
# R script to run FlowSOM and FlowSOM_meta
#
# Lukas M. Weber, December 2015
#########################################################################################


library(flowCore)
library(FlowSOM)



#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare.fcs")

# FlowSOM requires input data as flowFrame objects

data_Levine_32 <- flowCore::read.FCS(file_Levine_32, transformation = FALSE)
data_Levine_13 <- flowCore::read.FCS(file_Levine_13, transformation = FALSE)
data_Mosmann <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)

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



###################
### Run FlowSOM ###
###################

# run FlowSOM

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

set.seed(123)
runtime_FlowSOM_Mosmann <- system.time({
  fSOM_Mosmann <- FlowSOM::ReadInput(data_Mosmann, transform = FALSE, scale = FALSE)
  fSOM_Mosmann <- FlowSOM::BuildSOM(fSOM_Mosmann, colsToUse = marker_cols_Mosmann)
  fSOM_Mosmann <- FlowSOM::BuildMST(fSOM_Mosmann)
})

# plots

FlowSOM::PlotStars(fSOM_Levine_32)
FlowSOM::PlotStars(fSOM_Levine_13)
FlowSOM::PlotStars(fSOM_Mosmann)

# extract cluster labels

str(fSOM_Levine_32$map)

head(fSOM_Levine_32$map$mapping)
dim(fSOM_Levine_32$map$mapping)

clus_FlowSOM_Levine_32 <- fSOM_Levine_32$map$mapping[, 1]
clus_FlowSOM_Levine_13 <- fSOM_Levine_13$map$mapping[, 1]
clus_FlowSOM_Mosmann <- fSOM_Mosmann$map$mapping[, 1]

length(clus_FlowSOM_Levine_32)
length(clus_FlowSOM_Levine_13)
length(clus_FlowSOM_Mosmann)

# cluster sizes and number of clusters

table(clus_FlowSOM_Levine_32)
table(clus_FlowSOM_Levine_13)
table(clus_FlowSOM_Mosmann)

length(table(clus_FlowSOM_Levine_32))
length(table(clus_FlowSOM_Levine_13))
length(table(clus_FlowSOM_Mosmann))



########################
### Run FlowSOM_meta ###
########################

# run optional metaclustering step (FlowSOM_meta)

# set number of clusters
k_Levine_32 <- 20
k_Levine_13 <- 30
k_Mosmann <- 50

set.seed(123)
runtime_FlowSOM_meta_Levine_32 <- system.time({
  meta_clustering_Levine_32 <- FlowSOM::metaClustering_consensus(fSOM_Levine_32$map$codes, k = k_Levine_32)
})

set.seed(123)
runtime_FlowSOM_meta_Levine_13 <- system.time({
  meta_clustering_Levine_13 <- FlowSOM::metaClustering_consensus(fSOM_Levine_13$map$codes, k = k_Levine_13)
})

set.seed(123)
runtime_FlowSOM_meta_Mosmann <- system.time({
  meta_clustering_Mosmann <- FlowSOM::metaClustering_consensus(fSOM_Mosmann$map$codes, k = k_Mosmann)
})


# alternatively: set number of clusters automatically (does not perform well)

# set.seed(123)
# runtime_FlowSOM_meta_Levine_32 <- system.time({
# meta_clustering_Levine_32 <- FlowSOM::MetaClustering(fSOM_Levine_32$map$codes, 
#                                                      method = "metaClustering_consensus", 
#                                                      max = 50)
# })
# 
# set.seed(123)
# runtime_FlowSOM_meta_Levine_13 <- system.time({
# meta_clustering_Levine_13 <- FlowSOM::MetaClustering(fSOM_Levine_13$map$codes, 
#                                                      method = "metaClustering_consensus", 
#                                                      max = 50)
# })
# 
# set.seed(123)
# runtime_FlowSOM_meta_Mosmann <- system.time({
#   meta_clustering_Mosmann <- FlowSOM::MetaClustering(fSOM_Mosmann$map$codes, 
#                                                     method = "metaClustering_consensus", 
#                                                     max = 50)
# })


# combine runtime

runtime_FlowSOM_meta_Levine_32 <- runtime_FlowSOM_Levine_32 + runtime_FlowSOM_meta_Levine_32
runtime_FlowSOM_meta_Levine_13 <- runtime_FlowSOM_Levine_13 + runtime_FlowSOM_meta_Levine_13
runtime_FlowSOM_meta_Mosmann <- runtime_FlowSOM_Mosmann + runtime_FlowSOM_meta_Mosmann

# extract cluster labels

meta_clustering_Levine_32
meta_clustering_Levine_13
meta_clustering_Mosmann

clus_FlowSOM_meta_Levine_32 <- meta_clustering_Levine_32[fSOM_Levine_32$map$mapping[, 1]]
clus_FlowSOM_meta_Levine_13 <- meta_clustering_Levine_13[fSOM_Levine_13$map$mapping[, 1]]
clus_FlowSOM_meta_Mosmann <- meta_clustering_Mosmann[fSOM_Mosmann$map$mapping[, 1]]

length(clus_FlowSOM_meta_Levine_32)
length(clus_FlowSOM_meta_Levine_13)
length(clus_FlowSOM_meta_Mosmann)

# cluster sizes and number of clusters

table(clus_FlowSOM_meta_Levine_32)
table(clus_FlowSOM_meta_Levine_13)
table(clus_FlowSOM_meta_Mosmann)

length(table(clus_FlowSOM_meta_Levine_32))
length(table(clus_FlowSOM_meta_Levine_13))
length(table(clus_FlowSOM_meta_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_FlowSOM_Levine_32 <- data.frame(label = clus_FlowSOM_Levine_32)
res_FlowSOM_Levine_13 <- data.frame(label = clus_FlowSOM_Levine_13)
res_FlowSOM_Mosmann <- data.frame(label = clus_FlowSOM_Mosmann)

res_FlowSOM_meta_Levine_32 <- data.frame(label = clus_FlowSOM_meta_Levine_32)
res_FlowSOM_meta_Levine_13 <- data.frame(label = clus_FlowSOM_meta_Levine_13)
res_FlowSOM_meta_Mosmann <- data.frame(label = clus_FlowSOM_meta_Mosmann)


write.table(res_FlowSOM_Levine_32, 
            file = "../results/FlowSOM/FlowSOM_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Levine_13, 
            file = "../results/FlowSOM/FlowSOM_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_Mosmann, 
            file = "../results/FlowSOM/FlowSOM_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_FlowSOM_meta_Levine_32, 
            file = "../results/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Levine_13, 
            file = "../results/FlowSOM_meta/FlowSOM_meta_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_FlowSOM_meta_Mosmann, 
            file = "../results/FlowSOM_meta/FlowSOM_meta_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_FlowSOM <- t(data.frame(
  Levine_2015_marrow_32 = runtime_FlowSOM_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_FlowSOM_Levine_13["elapsed"], 
  Mosmann_2014_rare = runtime_FlowSOM_Mosmann["elapsed"], 
  row.names = "runtime"))

runtime_FlowSOM_meta <- t(data.frame(
  Levine_2015_marrow_32 = runtime_FlowSOM_meta_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_FlowSOM_meta_Levine_13["elapsed"], 
  Mosmann_2014_rare = runtime_FlowSOM_meta_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_FlowSOM, file = "../results/runtime/runtime_FlowSOM.txt", quote = FALSE, sep = "\t")

write.table(runtime_FlowSOM_meta, file = "../results/runtime/runtime_FlowSOM_meta.txt", quote = FALSE, sep = "\t")

# save session information

sink(file = "../results/session_info/FlowSOM_and_FlowSOM_meta_session_info.txt")
sessionInfo()
sink()

# save R objects

save.image(file = "../results/RData_files/FlowSOM_and_FlowSOM_meta_results.RData")


