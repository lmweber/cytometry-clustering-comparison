#########################################################################################
# R script to run k-means
#
# Lukas M. Weber, October 2015
#########################################################################################


library(flowCore)


# load data

file <- "../data/Levine_BMMC_32/Levine_BMMC_32.fcs"
data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE))

data_kmeans <- data[, -grep("label", colnames(data))]

head(data_kmeans)
dim(data_kmeans)


# run k-means

k <- 20  # number of clusters

set.seed(1234)

system.time(
out_kmeans <- kmeans(data_kmeans, k)
)


# extract cluster labels

clus_kmeans <- out_kmeans$cluster

table(clus_kmeans)
length(clus_kmeans)

res <- data.frame(label = clus_kmeans)
write.table(res, file = "../results/Levine_BMMC_32/kmeans/Levine_BMMC_32_kmeans_labels.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

