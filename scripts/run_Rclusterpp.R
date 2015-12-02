#########################################################################################
# R script to run Rclusterpp
#
# Lukas M. Weber, December 2015
#########################################################################################


library(flowCore)
library(Rclusterpp)


# note: Rclusterpp did not complete after 3 days with 16 cores for the Mosmann_2014_rare
# data set. Code for this data set has been commented out below.


#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13.fcs")
# file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare.fcs")

data_Levine_32 <- flowCore::exprs(flowCore::read.FCS(file_Levine_32, transformation = FALSE))
data_Levine_13 <- flowCore::exprs(flowCore::read.FCS(file_Levine_13, transformation = FALSE))
# data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

head(data_Levine_32)
head(data_Levine_13)
# head(data_Mosmann)

dim(data_Levine_32)
dim(data_Levine_13)
# dim(data_Mosmann)

# indices of protein marker columns

marker_cols_Levine_32 <- 5:36
marker_cols_Levine_13 <- 1:13
# marker_cols_Mosmann <- 7:21

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
# length(marker_cols_Mosmann)

# subset data

data_Levine_32 <- data_Levine_32[, marker_cols_Levine_32]
data_Levine_13 <- data_Levine_13[, marker_cols_Levine_13]
# data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine_32)
dim(data_Levine_13)
# dim(data_Mosmann)



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

# set.seed(123)
# Rclusterpp.setThreads(n_cores)  # set number of cores
# runtime_Mosmann <- system.time({
#   out_Rclusterpp_Mosmann <- Rclusterpp.hclust(data_Mosmann, method = "average", distance = "euclidean")
# })


# cut dendrogram at an arbitrary number of clusters k and extract cluster labels

k_Levine_32 <- 20
k_Levine_13 <- 30
# k_Mosmann <- 50

clus_Rclusterpp_Levine_32 <- cutree(out_Rclusterpp_Levine_32, k = k_Levine_32)
clus_Rclusterpp_Levine_13 <- cutree(out_Rclusterpp_Levine_13, k = k_Levine_13)
# clus_Rclusterpp_Mosmann <- cutree(out_Rclusterpp_Mosmann, k = k_Mosmann)

length(clus_Rclusterpp_Levine_32)
length(clus_Rclusterpp_Levine_13)
# length(clus_Rclusterpp_Mosmann)


# cluster sizes and number of clusters

table(clus_Rclusterpp_Levine_32)
table(clus_Rclusterpp_Levine_13)
# table(clus_Rclusterpp_Mosmann)

length(table(clus_Rclusterpp_Levine_32))
length(table(clus_Rclusterpp_Levine_13))
# length(table(clus_Rclusterpp_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_Rclusterpp_Levine_32 <- data.frame(label = clus_Rclusterpp_Levine_32)
res_Rclusterpp_Levine_13 <- data.frame(label = clus_Rclusterpp_Levine_13)
# res_Rclusterpp_Mosmann <- data.frame(label = clus_Rclusterpp_Mosmann)

write.table(res_Rclusterpp_Levine_32, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_Rclusterpp_Levine_13, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
# write.table(res_Rclusterpp_Mosmann, 
#             file = "../results/Rclusterpp/Rclusterpp_labels_Mosmann_2014_rare.txt", 
#             row.names = FALSE, quote = FALSE, sep = "\t")

# save runtime

runtime_Rclusterpp <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_Levine_13["elapsed"], 
#   Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_Rclusterpp, file = "../results/runtime/runtime_Rclusterpp.txt", quote = FALSE, sep = "\t")

# save session information

sink(file = "../results/session_info/Rclusterpp_session_info.txt")
sessionInfo()
sink()

# save R objects

save.image(file = "../results/RData_files/Rclusterpp_results.RData")


