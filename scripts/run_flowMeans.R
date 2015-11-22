#########################################################################################
# R script to run flowMeans
#
# Lukas M. Weber, November 2015
#########################################################################################


library(flowCore)
library(flowMeans)



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



#####################
### Run flowMeans ###
#####################

# run flowMeans

set.seed(123)
runtime_Levine <- system.time({
  out_flowMeans_Levine <- flowMeans(data_Levine, Standardize = FALSE)
})

set.seed(123)
runtime_Mosmann <- system.time({
  out_flowMeans_Mosmann <- flowMeans(data_Mosmann, Standardize = FALSE, NumC = 50)
})

# extract cluster labels

clus_flowMeans_Levine <- out_flowMeans_Levine@Label
clus_flowMeans_Mosmann <- out_flowMeans_Mosmann@Label

length(clus_flowMeans_Levine)
length(clus_flowMeans_Mosmann)

# cluster sizes and number of clusters

table(clus_flowMeans_Levine)
table(clus_flowMeans_Mosmann)

length(table(clus_flowMeans_Levine))
length(table(clus_flowMeans_Mosmann))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_flowMeans_Levine <- data.frame(label = clus_flowMeans_Levine)
res_flowMeans_Mosmann <- data.frame(label = clus_flowMeans_Mosmann)

write.table(res_flowMeans_Levine, 
            file = "../results/flowMeans/flowMeans_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_flowMeans_Mosmann, 
            file = "../results/flowMeans/flowMeans_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# save runtime

runtime_flowMeans <- t(data.frame(
  Levine_2015_marrow_32 = runtime_Levine["elapsed"], 
  Mosmann_2014_rare = runtime_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_flowMeans, file = "../results/runtime/runtime_flowMeans.txt", quote = FALSE, sep = "\t")

# save session information

sink(file = "../results/session_info/flowMeans_session_info.txt")
sessionInfo()
sink()

# save R objects

save.image(file = "../results/RData_files/flowMeans_results.RData")


