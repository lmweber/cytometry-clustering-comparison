#########################################################################################
# R script to run Rclusterpp
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)
library(Rclusterpp)


# note: Rclusterpp did not complete after 3 days with 16 cores for the Mosmann_2014_rare
# data set. Code for this data set has been commented out below.


#################
### LOAD DATA ###
#################

DATA_DIR <- "../../benchmark_data_sets"

file_Levine <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32.fcs")
# file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare.fcs")

data_Levine <- flowCore::exprs(flowCore::read.FCS(file_Levine, transformation = FALSE))
# data_Mosmann <- flowCore::exprs(flowCore::read.FCS(file_Mosmann, transformation = FALSE))

head(data_Levine)
# head(data_Mosmann)

dim(data_Levine)
# dim(data_Mosmann)

# indices of protein marker columns

marker_cols_Levine <- 5:36
# marker_cols_Mosmann <- 7:21

length(marker_cols_Levine)
# length(marker_cols_Mosmann)

# subset data

data_Levine <- data_Levine[, marker_cols_Levine]
# data_Mosmann <- data_Mosmann[, marker_cols_Mosmann]

dim(data_Levine)
# dim(data_Mosmann)



######################
### Run Rclusterpp ###
######################

# run Rclusterpp
# note: uses maximum number of cores if setThreads is left as default

n_cores <- 16

Rclusterpp.setThreads(n_cores)  # set number of cores
set.seed(123)
runtime_Levine <- system.time({
  out_Rclusterpp_Levine <- Rclusterpp.hclust(data_Levine, method = "average", distance = "euclidean")
})

# Rclusterpp.setThreads(n_cores)  # set number of cores
# set.seed(123)
# runtime_Mosmann <- system.time({
#   out_Rclusterpp_Mosmann <- Rclusterpp.hclust(data_Mosmann, method = "average", distance = "euclidean")
# })


# cut dendrogram at an arbitrary number of clusters k and extract cluster labels

k_Levine <- 20
# k_Mosmann <- 50

clus_Rclusterpp_Levine <- cutree(out_Rclusterpp_Levine, k = k_Levine)
# clus_Rclusterpp_Mosmann <- cutree(out_Rclusterpp_Mosmann, k = k_Mosmann)

length(clus_Rclusterpp_Levine)
# length(clus_Rclusterpp_Mosmann)


# cluster sizes and number of clusters

table(clus_Rclusterpp_Levine)
# table(clus_Rclusterpp_Mosmann)

length(table(clus_Rclusterpp_Levine))
# length(table(clus_Rclusterpp_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_Rclusterpp_Levine <- data.frame(label = clus_Rclusterpp_Levine)
# res_Rclusterpp_Mosmann <- data.frame(label = clus_Rclusterpp_Mosmann)

write.table(res_Rclusterpp_Levine, 
            file = "../results/Rclusterpp/Rclusterpp_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
# write.table(res_Rclusterpp_Mosmann, 
#             file = "../results/Rclusterpp/Rclusterpp_labels_Mosmann_2014_rare.txt", 
#             row.names = FALSE, quote = FALSE, sep = "\t")

# save runtime

runtime_Rclusterpp <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine["elapsed"], 
#   Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_Rclusterpp, file = "../results/runtime/runtime_Rclusterpp.txt", quote = FALSE, sep = "\t")

# save session information

sink(file = "../results/session_info/Rclusterpp_session_info.txt")
sessionInfo()
sink()

# save R objects

save.image(file = "../results/RData_files/Rclusterpp_results.RData")


