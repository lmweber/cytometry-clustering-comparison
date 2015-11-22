#########################################################################################
# R script to run k-means
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)



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



###################
### Run k-means ###
###################

# number of clusters

k_Levine <- 20
k_Mosmann <- 50


# run k-means
# note: does not converge for some random seeds

set.seed(1234)
runtime_Levine <- system.time({
  out_kmeans_Levine <- kmeans(data_Levine, k_Levine)
})

set.seed(1000)
runtime_Mosmann <- system.time({
  out_kmeans_Mosmann <- kmeans(data_Mosmann, k_Mosmann)
})


# extract cluster labels

clus_kmeans_Levine <- out_kmeans_Levine$cluster
clus_kmeans_Mosmann <- out_kmeans_Mosmann$cluster

length(clus_kmeans_Levine)
length(clus_kmeans_Mosmann)


# cluster sizes and number of clusters

table(clus_kmeans_Levine)
table(clus_kmeans_Mosmann)

length(table(clus_kmeans_Levine))
length(table(clus_kmeans_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_kmeans_Levine <- data.frame(label = clus_kmeans_Levine)
res_kmeans_Mosmann <- data.frame(label = clus_kmeans_Mosmann)

write.table(res_kmeans_Levine, 
            file = "../results/kmeans/kmeans_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_kmeans_Mosmann, 
            file = "../results/kmeans/kmeans_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_kmeans <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine["elapsed"], 
  Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_kmeans, file = "../results/runtime/runtime_kmeans.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/kmeans_session_info.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/kmeans_results.RData")


