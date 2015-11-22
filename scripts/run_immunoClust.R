#########################################################################################
# R script to run immunoClust
#
# Lukas M. Weber, November 2015
#########################################################################################


# note installation from Bioconductor requires GNU Scientific Library

library(flowCore)
library(immunoClust)



#################
### LOAD DATA ###
#################

# use non-transformed data files, since immunoClust will transform automatically

DATA_DIR <- "../../benchmark_data_sets"

file_Levine <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransf.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_rare/data/Mosmann_2014_rare_notransf.fcs")

# input data as flowFrame objects

data_Levine <- flowCore::read.FCS(file_Levine, transformation = FALSE)
data_Mosmann <- flowCore::read.FCS(file_Mosmann, transformation = FALSE)

head(data_Levine)
head(data_Mosmann)

dim(data_Levine)
dim(data_Mosmann)

# indices of protein marker columns

marker_cols_Levine <- 5:36
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine)
length(marker_cols_Mosmann)

# column names (parameters)

pars_Levine <- colnames(data_Levine)[marker_cols_Levine]
pars_Mosmann <- colnames(data_Mosmann)[marker_cols_Mosmann]

length(pars_Levine)
length(pars_Mosmann)



###########################################
### Run immunoClust and immunoClust_all ###
###########################################

# run immunoClust
# (note: decreasing the bias argument increases the number of clusters)

runtime_immunoClust_Levine <- system.time({
  out_immunoClust_Levine <- immunoClust::cell.process(data_Levine, parameters = pars_Levine)
})

runtime_immunoClust_Mosmann <- system.time({
  out_immunoClust_Mosmann <- immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, 
                                                       bias = 0.1)
})


# run immmunoClust_all (with additional step to classify all cells)

runtime_immunoClust_all_Levine <- system.time({
  out_immunoClust_all_Levine <- immunoClust::cell.process(data_Levine, parameters = pars_Levine, 
                                                          classify.all = TRUE)
})

runtime_immunoClust_all_Mosmann <- system.time({
  out_immunoClust_all_Mosmann <- immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, 
                                                           bias = 0.1, classify.all = TRUE)
})


# extract cluster labels

# immunoClust
summary(out_immunoClust_Levine)  # number of clusters
summary(out_immunoClust_Mosmann)

clus_immunoClust_Levine <- out_immunoClust_Levine@label
clus_immunoClust_Mosmann <- out_immunoClust_Mosmann@label

# immunoClust_all
summary(out_immunoClust_all_Levine)
summary(out_immunoClust_all_Mosmann)

clus_immunoClust_all_Levine <- out_immunoClust_all_Levine@label
clus_immunoClust_all_Mosmann <- out_immunoClust_all_Mosmann@label


# cluster sizes and number of clusters

# immunoClust
table(clus_immunoClust_Levine)
table(clus_immunoClust_Mosmann)

length(table(clus_immunoClust_Levine))
length(table(clus_immunoClust_Mosmann))

# immunoClust_all
table(clus_immunoClust_all_Levine)
table(clus_immunoClust_all_Mosmann)

length(table(clus_immunoClust_all_Levine))
length(table(clus_immunoClust_all_Mosmann))


# plots

# immunoClust
data_transf_Levine <- immunoClust::trans.ApplyToData(out_immunoClust_Levine, data_Levine)
data_transf_Mosmann <- immunoClust::trans.ApplyToData(out_immunoClust_Mosmann, data_Mosmann)

png("../results/immunoClust/plot_immunoClust_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine, data_transf_Levine, N = 1000)
dev.off()

png("../results/immunoClust/plot_immunoClust_Mosmann_2014_rare.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Mosmann, data_transf_Mosmann, N = 1000)
dev.off()

# immunoClust_all
data_transf_all_Levine <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine, data_Levine)
data_transf_all_Mosmann <- immunoClust::trans.ApplyToData(out_immunoClust_all_Mosmann, data_Mosmann)

png("../results/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine, data_transf_all_Levine, N = 1000)
dev.off()

png("../results/immunoClust_all/plot_immunoClust_all_Mosmann_2014_rare.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Mosmann, data_transf_all_Mosmann, N = 1000)
dev.off()



####################
### SAVE RESULTS ###
####################

# save cluster labels

res_immunoClust_Levine <- data.frame(label = clus_immunoClust_Levine)
res_immunoClust_Mosmann <- data.frame(label = clus_immunoClust_Mosmann)

res_immunoClust_all_Levine <- data.frame(label = clus_immunoClust_all_Levine)
res_immunoClust_all_Mosmann <- data.frame(label = clus_immunoClust_all_Mosmann)


write.table(res_immunoClust_Levine, 
            file = "../results/immunoClust/immunoClust_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Mosmann, 
            file = "../results/immunoClust/immunoClust_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_immunoClust_all_Levine, 
            file = "../results/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Mosmann, 
            file = "../results/immunoClust_all/immunoClust_all_labels_Mosmann_2014_rare.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_immunoClust <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_Levine["elapsed"], 
  Mosmann_2014_rare = runtime_immunoClust_Mosmann["elapsed"], 
  row.names = "runtime"))

runtime_immunoClust_all <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_all_Levine["elapsed"], 
  Mosmann_2014_rare = runtime_immunoClust_all_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_immunoClust, file = "../results/runtime/runtime_immunoClust.txt", quote = FALSE, sep = "\t")

write.table(runtime_immunoClust_all, file = "../results/runtime/runtime_immunoClust_all.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/immunoClust_and_immunoClust_all_session_info.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/immunoClust_and_immunoClust_all_results.RData")


