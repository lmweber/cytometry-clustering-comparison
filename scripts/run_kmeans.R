#########################################################################################
# R script to run k-means
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)



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




#########################################################
### Run k-means: manually selected number of clusters ###
#########################################################

# number of clusters
k_Levine_32 <- 40
k_Levine_13 <- 40
k_Nilsson <- 40
k_Mosmann <- 40


# run k-means
# note: returns errors for some random seeds; also additional iterations required

set.seed(1234)
runtime_Levine_32 <- system.time({
  out_kmeans_Levine_32 <- kmeans(data_Levine_32, k_Levine_32)
})

set.seed(1000)
runtime_Levine_13 <- system.time({
  out_kmeans_Levine_13 <- kmeans(data_Levine_13, k_Levine_13)
})

set.seed(1234)
runtime_Nilsson <- system.time({
  out_kmeans_Nilsson <- kmeans(data_Nilsson, k_Nilsson, iter.max = 50)
})

set.seed(2000)
runtime_Mosmann <- system.time({
  out_kmeans_Mosmann <- kmeans(data_Mosmann, k_Mosmann, iter.max = 50)
})


# extract cluster labels

clus_kmeans_Levine_32 <- out_kmeans_Levine_32$cluster
clus_kmeans_Levine_13 <- out_kmeans_Levine_13$cluster
clus_kmeans_Nilsson <- out_kmeans_Nilsson$cluster
clus_kmeans_Mosmann <- out_kmeans_Mosmann$cluster

length(clus_kmeans_Levine_32)
length(clus_kmeans_Levine_13)
length(clus_kmeans_Nilsson)
length(clus_kmeans_Mosmann)


# cluster sizes and number of clusters

table(clus_kmeans_Levine_32)
table(clus_kmeans_Levine_13)
table(clus_kmeans_Nilsson)
table(clus_kmeans_Mosmann)

length(table(clus_kmeans_Levine_32))
length(table(clus_kmeans_Levine_13))
length(table(clus_kmeans_Nilsson))
length(table(clus_kmeans_Mosmann))


# save cluster labels

res_kmeans_Levine_32 <- data.frame(label = clus_kmeans_Levine_32)
res_kmeans_Levine_13 <- data.frame(label = clus_kmeans_Levine_13)
res_kmeans_Nilsson <- data.frame(label = clus_kmeans_Nilsson)
res_kmeans_Mosmann <- data.frame(label = clus_kmeans_Mosmann)

write.table(res_kmeans_Levine_32, 
            file = "../results_manual/kmeans/kmeans_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_kmeans_Levine_13, 
            file = "../results_manual/kmeans/kmeans_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_kmeans_Nilsson, 
            file = "../results_manual/kmeans/kmeans_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_kmeans_Mosmann, 
            file = "../results_manual/kmeans/kmeans_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_kmeans <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_kmeans, file = "../results_manual/runtime/runtime_kmeans.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_manual/session_info/session_info_kmeans.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_manual/RData_files/results_kmeans.RData")


