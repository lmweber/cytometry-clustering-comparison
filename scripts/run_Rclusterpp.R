#########################################################################################
# R script to run Rclusterpp
#
# Lukas M. Weber, March 2016
#########################################################################################


library(flowCore)
library(Rclusterpp)


# note: subsampling is required for larger data sets due to slow runtime



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


# subsampling (required for larger data sets due to slow runtime)

n_sub_Levine_32 <- min(100000, nrow(data_Levine_32))
n_sub_Levine_13 <- min(100000, nrow(data_Levine_13))
n_sub_Nilsson <- min(100000, nrow(data_Nilsson))
n_sub_Mosmann <- min(100000, nrow(data_Mosmann))

data_Levine_32 <- data_Levine_32[sample(1:nrow(data_Levine_32), n_sub_Levine_32, replace = FALSE), ]
data_Levine_13 <- data_Levine_13[sample(1:nrow(data_Levine_13), n_sub_Levine_13, replace = FALSE), ]
data_Nilsson <- data_Nilsson[sample(1:nrow(data_Nilsson), n_sub_Nilsson, replace = FALSE), ]
data_Mosmann <- data_Mosmann[sample(1:nrow(data_Mosmann), n_sub_Mosmann, replace = FALSE), ]

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)


# save subsampled data (contains true population labels for subsampled rows)

flowCore::write.FCS(flowCore::flowFrame(data_Levine_32), filename = "../results/Rclusterpp/Levine_2015_marrow_32_sub.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_Levine_13), filename = "../results/Rclusterpp/Levine_2015_marrow_13_sub.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_Nilsson), filename = "../results/Rclusterpp/Nilsson_2013_HSC_sub.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_Mosmann), filename = "../results/Rclusterpp/Mosmann_2014_activ_sub.fcs")


# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
marker_cols_Nilsson <- c(5:7, 9:18)
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)


# subset columns

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
data_Nilsson <- data_Nilsson[, marker_cols_Nilsson]
data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine_32)
dim(data_Levine_13)
dim(data_Nilsson)
dim(data_Mosmann)




######################
### Run Rclusterpp ###
######################

# run Rclusterpp
# note: uses maximum number of cores if setThreads is left as default

n_cores <- 16


set.seed(123)
Rclusterpp.setThreads(n_cores)  # set number of cores
runtime_Levine_32 <- system.time({
  out_Rclusterpp_Levine_32 <- Rclusterpp.hclust(data_Levine_32, method = "average", distance = "euclidean")
})


set.seed(123)
Rclusterpp.setThreads(n_cores)  # set number of cores
runtime_Levine_13 <- system.time({
  out_Rclusterpp_Levine_13 <- Rclusterpp.hclust(data_Levine_13, method = "average", distance = "euclidean")
})


set.seed(123)
Rclusterpp.setThreads(n_cores)  # set number of cores
runtime_Nilsson <- system.time({
  out_Rclusterpp_Nilsson <- Rclusterpp.hclust(data_Nilsson, method = "average", distance = "euclidean")
})


set.seed(123)
Rclusterpp.setThreads(n_cores)  # set number of cores
runtime_Mosmann <- system.time({
  out_Rclusterpp_Mosmann <- Rclusterpp.hclust(data_Mosmann, method = "average", distance = "euclidean")
})


# cut dendrogram at an arbitrary number of clusters k and extract cluster labels

k_Levine_32 <- 40
k_Levine_13 <- 40
k_Nilsson <- 40
k_Mosmann <- 40

clus_Rclusterpp_Levine_32 <- cutree(out_Rclusterpp_Levine_32, k = k_Levine_32)
clus_Rclusterpp_Levine_13 <- cutree(out_Rclusterpp_Levine_13, k = k_Levine_13)
clus_Rclusterpp_Nilsson <- cutree(out_Rclusterpp_Nilsson, k = k_Nilsson)
clus_Rclusterpp_Mosmann <- cutree(out_Rclusterpp_Mosmann, k = k_Mosmann)

length(clus_Rclusterpp_Levine_32)
length(clus_Rclusterpp_Levine_13)
length(clus_Rclusterpp_Nilsson)
length(clus_Rclusterpp_Mosmann)


# cluster sizes and number of clusters

table(clus_Rclusterpp_Levine_32)
table(clus_Rclusterpp_Levine_13)
table(clus_Rclusterpp_Nilsson)
table(clus_Rclusterpp_Mosmann)

length(table(clus_Rclusterpp_Levine_32))
length(table(clus_Rclusterpp_Levine_13))
length(table(clus_Rclusterpp_Nilsson))
length(table(clus_Rclusterpp_Mosmann))




####################
### SAVE RESULTS ###
####################

# save cluster labels

res_Rclusterpp_Levine_32 <- data.frame(label = clus_Rclusterpp_Levine_32)
res_Rclusterpp_Levine_13 <- data.frame(label = clus_Rclusterpp_Levine_13)
res_Rclusterpp_Nilsson <- data.frame(label = clus_Rclusterpp_Nilsson)
res_Rclusterpp_Mosmann <- data.frame(label = clus_Rclusterpp_Mosmann)

write.table(res_Rclusterpp_Levine_32, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_Rclusterpp_Levine_13, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_Rclusterpp_Nilsson, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_Rclusterpp_Mosmann, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_Rclusterpp <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_Rclusterpp, file = "../results/runtime/runtime_Rclusterpp.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/session_info_Rclusterpp.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/results_Rclusterpp.RData")


