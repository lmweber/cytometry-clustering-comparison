#########################################################################################
# R script to run immunoClust and immunoClust_all
#
# Lukas M. Weber, March 2016
#########################################################################################


# note installation from Bioconductor requires GNU Scientific Library

library(flowCore)
library(immunoClust)



#################
### LOAD DATA ###
#################

# use non-transformed data files, since immunoClust will transform automatically

DATA_DIR <- "../../benchmark_data_sets"

file_Levine_32 <- file.path(DATA_DIR, "Levine_2015_marrow_32/data/Levine_2015_marrow_32_notransform.fcs")
file_Levine_13 <- file.path(DATA_DIR, "Levine_2015_marrow_13/data/Levine_2015_marrow_13_notransform.fcs")
file_Nilsson <- file.path(DATA_DIR, "Nilsson_2013_HSC/data/Nilsson_2013_HSC_notransform.fcs")
file_Mosmann <- file.path(DATA_DIR, "Mosmann_2014_activ/data/Mosmann_2014_activ_notransform.fcs")


# input data as flowFrame objects

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
marker_cols_Mosmann <- 7:21

length(marker_cols_Levine_32)
length(marker_cols_Levine_13)
length(marker_cols_Nilsson)
length(marker_cols_Mosmann)


# column names (parameters)

pars_Levine_32 <- colnames(data_Levine_32)[marker_cols_Levine_32]
pars_Levine_13 <- colnames(data_Levine_13)[marker_cols_Levine_13]
pars_Nilsson <- colnames(data_Nilsson)[marker_cols_Nilsson]
pars_Mosmann <- colnames(data_Mosmann)[marker_cols_Mosmann]

length(pars_Levine_32)
length(pars_Levine_13)
length(pars_Nilsson)
length(pars_Mosmann)




###########################################
### Run immunoClust and immunoClust_all ###
###########################################

# run immunoClust
# (note: decreasing the bias argument increases the number of clusters)

set.seed(123)
runtime_immunoClust_Levine_32 <- system.time({
  out_immunoClust_Levine_32 <- immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32)
})

set.seed(123)
runtime_immunoClust_Levine_13 <- system.time({
  out_immunoClust_Levine_13 <- immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13)
})

set.seed(123)
runtime_immunoClust_Nilsson <- system.time({
  out_immunoClust_Nilsson <- immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, bias = 0.1)
})

set.seed(123)
runtime_immunoClust_Mosmann <- system.time({
  out_immunoClust_Mosmann <- immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, bias = 0.1)
})


# run immmunoClust_all (with additional step to classify all cells)

set.seed(123)
runtime_immunoClust_all_Levine_32 <- system.time({
  out_immunoClust_all_Levine_32 <- immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32, 
                                                             classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Levine_13 <- system.time({
  out_immunoClust_all_Levine_13 <- immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13, 
                                                             classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Nilsson <- system.time({
  out_immunoClust_all_Nilsson <- immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, 
                                                           bias = 0.1, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Mosmann <- system.time({
  out_immunoClust_all_Mosmann <- immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, 
                                                           bias = 0.1, classify.all = TRUE)
})


# extract cluster labels

# immunoClust

summary(out_immunoClust_Levine_32)  # number of clusters
summary(out_immunoClust_Levine_13)
summary(out_immunoClust_Nilsson)
summary(out_immunoClust_Mosmann)

clus_immunoClust_Levine_32 <- out_immunoClust_Levine_32@label
clus_immunoClust_Levine_13 <- out_immunoClust_Levine_13@label
clus_immunoClust_Nilsson <- out_immunoClust_Nilsson@label
clus_immunoClust_Mosmann <- out_immunoClust_Mosmann@label


# immunoClust_all

summary(out_immunoClust_all_Levine_32)
summary(out_immunoClust_all_Levine_13)
summary(out_immunoClust_all_Nilsson)
summary(out_immunoClust_all_Mosmann)

clus_immunoClust_all_Levine_32 <- out_immunoClust_all_Levine_32@label
clus_immunoClust_all_Levine_13 <- out_immunoClust_all_Levine_13@label
clus_immunoClust_all_Nilsson <- out_immunoClust_all_Nilsson@label
clus_immunoClust_all_Mosmann <- out_immunoClust_all_Mosmann@label


# cluster sizes and number of clusters

# immunoClust

table(clus_immunoClust_Levine_32)
table(clus_immunoClust_Levine_13)
table(clus_immunoClust_Nilsson)
table(clus_immunoClust_Mosmann)

length(table(clus_immunoClust_Levine_32))
length(table(clus_immunoClust_Levine_13))
length(table(clus_immunoClust_Nilsson))
length(table(clus_immunoClust_Mosmann))

# immunoClust_all

table(clus_immunoClust_all_Levine_32)
table(clus_immunoClust_all_Levine_13)
table(clus_immunoClust_all_Nilsson)
table(clus_immunoClust_all_Mosmann)

length(table(clus_immunoClust_all_Levine_32))
length(table(clus_immunoClust_all_Levine_13))
length(table(clus_immunoClust_all_Nilsson))
length(table(clus_immunoClust_all_Mosmann))


# plots

# immunoClust

data_transf_Levine_32 <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_32, data_Levine_32)
data_transf_Levine_13 <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_13, data_Levine_13)
data_transf_Nilsson <- immunoClust::trans.ApplyToData(out_immunoClust_Nilsson, data_Nilsson)
data_transf_Mosmann <- immunoClust::trans.ApplyToData(out_immunoClust_Mosmann, data_Mosmann)

png("../results/immunoClust/plot_immunoClust_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine_32, data_transf_Levine_32, N = 1000)
dev.off()

png("../results/immunoClust/plot_immunoClust_Levine_2015_marrow_13.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine_13, data_transf_Levine_13, N = 1000)
dev.off()

png("../results/immunoClust/plot_immunoClust_Nilsson_2013_HSC.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Nilsson, data_transf_Nilsson, N = 1000)
dev.off()

png("../results/immunoClust/plot_immunoClust_Mosmann_2014_activ.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Mosmann, data_transf_Mosmann, N = 1000)
dev.off()


# immunoClust_all

data_transf_all_Levine_32 <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_32, data_Levine_32)
data_transf_all_Levine_13 <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_13, data_Levine_13)
data_transf_all_Nilsson <- immunoClust::trans.ApplyToData(out_immunoClust_all_Nilsson, data_Nilsson)
data_transf_all_Mosmann <- immunoClust::trans.ApplyToData(out_immunoClust_all_Mosmann, data_Mosmann)

png("../results/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine_32, data_transf_all_Levine_32, N = 1000)
dev.off()

png("../results/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_13.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine_13, data_transf_all_Levine_13, N = 1000)
dev.off()

png("../results/immunoClust_all/plot_immunoClust_all_Nilsson_2013_HSC.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Nilsson, data_transf_all_Nilsson, N = 1000)
dev.off()

png("../results/immunoClust_all/plot_immunoClust_all_Mosmann_2014_activ.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Mosmann, data_transf_all_Mosmann, N = 1000)
dev.off()




####################
### SAVE RESULTS ###
####################

# save cluster labels

res_immunoClust_Levine_32 <- data.frame(label = clus_immunoClust_Levine_32)
res_immunoClust_Levine_13 <- data.frame(label = clus_immunoClust_Levine_13)
res_immunoClust_Nilsson <- data.frame(label = clus_immunoClust_Nilsson)
res_immunoClust_Mosmann <- data.frame(label = clus_immunoClust_Mosmann)

res_immunoClust_all_Levine_32 <- data.frame(label = clus_immunoClust_all_Levine_32)
res_immunoClust_all_Levine_13 <- data.frame(label = clus_immunoClust_all_Levine_13)
res_immunoClust_all_Nilsson <- data.frame(label = clus_immunoClust_all_Nilsson)
res_immunoClust_all_Mosmann <- data.frame(label = clus_immunoClust_all_Mosmann)


write.table(res_immunoClust_Levine_32, 
            file = "../results/immunoClust/immunoClust_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Levine_13, 
            file = "../results/immunoClust/immunoClust_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Nilsson, 
            file = "../results/immunoClust/immunoClust_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Mosmann, 
            file = "../results/immunoClust/immunoClust_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_immunoClust_all_Levine_32, 
            file = "../results/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Levine_13, 
            file = "../results/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Nilsson, 
            file = "../results/immunoClust_all/immunoClust_all_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Mosmann, 
            file = "../results/immunoClust_all/immunoClust_all_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_immunoClust <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_immunoClust_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_immunoClust_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_immunoClust_Mosmann["elapsed"], 
  row.names = "runtime"))

runtime_immunoClust_all <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_all_Levine_32["elapsed"], 
  Levine_2015_marrow_13 = runtime_immunoClust_all_Levine_13["elapsed"], 
  Nilsson_2013_HSC = runtime_immunoClust_all_Nilsson["elapsed"], 
  Mosmann_2014_activ = runtime_immunoClust_all_Mosmann["elapsed"], 
  row.names = "runtime"))

write.table(runtime_immunoClust, file = "../results/runtime/runtime_immunoClust.txt", quote = FALSE, sep = "\t")

write.table(runtime_immunoClust_all, file = "../results/runtime/runtime_immunoClust_all.txt", quote = FALSE, sep = "\t")


# save session information

sink(file = "../results/session_info/session_info_immunoClust_and_immunoClust_all.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results/RData_files/results_immunoClust_and_immunoClust_all.RData")


