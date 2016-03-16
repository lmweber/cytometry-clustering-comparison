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




#########################################################################
### Run immunoClust and immunoClust_all: automatic number of clusters ###
#########################################################################

# run immunoClust with automatic selection of number of clusters
# (note: decreasing the bias argument increases the number of clusters)

set.seed(123)
runtime_immunoClust_Levine_32_auto <- system.time({
  out_immunoClust_Levine_32_auto <- 
    immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32)
})

set.seed(123)
runtime_immunoClust_Levine_13_auto <- system.time({
  out_immunoClust_Levine_13_auto <- 
    immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13)
})

set.seed(123)
runtime_immunoClust_Nilsson_auto <- system.time({
  out_immunoClust_Nilsson_auto <- 
    immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson)
})

set.seed(123)
runtime_immunoClust_Mosmann_auto <- system.time({
  out_immunoClust_Mosmann_auto <- 
    immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann)
})


# run immmunoClust_all (additional step to classify all cells)

set.seed(123)
runtime_immunoClust_all_Levine_32_auto <- system.time({
  out_immunoClust_all_Levine_32_auto <- 
    immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Levine_13_auto <- system.time({
  out_immunoClust_all_Levine_13_auto <- 
    immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Nilsson_auto <- system.time({
  out_immunoClust_all_Nilsson_auto <- 
    immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Mosmann_auto <- system.time({
  out_immunoClust_all_Mosmann_auto <- 
    immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, classify.all = TRUE)
})


# extract cluster labels

# immunoClust

summary(out_immunoClust_Levine_32_auto)  # number of clusters
summary(out_immunoClust_Levine_13_auto)
summary(out_immunoClust_Nilsson_auto)
summary(out_immunoClust_Mosmann_auto)

clus_immunoClust_Levine_32_auto <- out_immunoClust_Levine_32_auto@label
clus_immunoClust_Levine_13_auto <- out_immunoClust_Levine_13_auto@label
clus_immunoClust_Nilsson_auto <- out_immunoClust_Nilsson_auto@label
clus_immunoClust_Mosmann_auto <- out_immunoClust_Mosmann_auto@label


# immunoClust_all

summary(out_immunoClust_all_Levine_32_auto)
summary(out_immunoClust_all_Levine_13_auto)
summary(out_immunoClust_all_Nilsson_auto)
summary(out_immunoClust_all_Mosmann_auto)

clus_immunoClust_all_Levine_32_auto <- out_immunoClust_all_Levine_32_auto@label
clus_immunoClust_all_Levine_13_auto <- out_immunoClust_all_Levine_13_auto@label
clus_immunoClust_all_Nilsson_auto <- out_immunoClust_all_Nilsson_auto@label
clus_immunoClust_all_Mosmann_auto <- out_immunoClust_all_Mosmann_auto@label


# cluster sizes and number of clusters

# immunoClust

table(clus_immunoClust_Levine_32_auto)
table(clus_immunoClust_Levine_13_auto)
table(clus_immunoClust_Nilsson_auto)
table(clus_immunoClust_Mosmann_auto)

length(table(clus_immunoClust_Levine_32_auto))
length(table(clus_immunoClust_Levine_13_auto))
length(table(clus_immunoClust_Nilsson_auto))
length(table(clus_immunoClust_Mosmann_auto))

# immunoClust_all

table(clus_immunoClust_all_Levine_32_auto)
table(clus_immunoClust_all_Levine_13_auto)
table(clus_immunoClust_all_Nilsson_auto)
table(clus_immunoClust_all_Mosmann_auto)

length(table(clus_immunoClust_all_Levine_32_auto))
length(table(clus_immunoClust_all_Levine_13_auto))
length(table(clus_immunoClust_all_Nilsson_auto))
length(table(clus_immunoClust_all_Mosmann_auto))


# plots

# immunoClust

data_transf_Levine_32_auto <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_32_auto, data_Levine_32)
data_transf_Levine_13_auto <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_13_auto, data_Levine_13)
data_transf_Nilsson_auto <- immunoClust::trans.ApplyToData(out_immunoClust_Nilsson_auto, data_Nilsson)
data_transf_Mosmann_auto <- immunoClust::trans.ApplyToData(out_immunoClust_Mosmann_auto, data_Mosmann)

png("../results_auto/immunoClust/plot_immunoClust_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine_32_auto, data_transf_Levine_32_auto, N = 1000)
dev.off()

png("../results_auto/immunoClust/plot_immunoClust_Levine_2015_marrow_13.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine_13_auto, data_transf_Levine_13_auto, N = 1000)
dev.off()

png("../results_auto/immunoClust/plot_immunoClust_Nilsson_2013_HSC.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Nilsson_auto, data_transf_Nilsson_auto, N = 1000)
dev.off()

png("../results_auto/immunoClust/plot_immunoClust_Mosmann_2014_activ.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Mosmann_auto, data_transf_Mosmann_auto, N = 1000)
dev.off()


# immunoClust_all

data_transf_all_Levine_32_auto <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_32_auto, data_Levine_32)
data_transf_all_Levine_13_auto <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_13_auto, data_Levine_13)
data_transf_all_Nilsson_auto <- immunoClust::trans.ApplyToData(out_immunoClust_all_Nilsson_auto, data_Nilsson)
data_transf_all_Mosmann_auto <- immunoClust::trans.ApplyToData(out_immunoClust_all_Mosmann_auto, data_Mosmann)

png("../results_auto/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine_32_auto, data_transf_all_Levine_32_auto, N = 1000)
dev.off()

png("../results_auto/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_13.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine_13_auto, data_transf_all_Levine_13_auto, N = 1000)
dev.off()

png("../results_auto/immunoClust_all/plot_immunoClust_all_Nilsson_2013_HSC.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Nilsson_auto, data_transf_all_Nilsson_auto, N = 1000)
dev.off()

png("../results_auto/immunoClust_all/plot_immunoClust_all_Mosmann_2014_activ.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Mosmann_auto, data_transf_all_Mosmann_auto, N = 1000)
dev.off()


# save cluster labels

res_immunoClust_Levine_32_auto <- data.frame(label = clus_immunoClust_Levine_32_auto)
res_immunoClust_Levine_13_auto <- data.frame(label = clus_immunoClust_Levine_13_auto)
res_immunoClust_Nilsson_auto <- data.frame(label = clus_immunoClust_Nilsson_auto)
res_immunoClust_Mosmann_auto <- data.frame(label = clus_immunoClust_Mosmann_auto)

res_immunoClust_all_Levine_32_auto <- data.frame(label = clus_immunoClust_all_Levine_32_auto)
res_immunoClust_all_Levine_13_auto <- data.frame(label = clus_immunoClust_all_Levine_13_auto)
res_immunoClust_all_Nilsson_auto <- data.frame(label = clus_immunoClust_all_Nilsson_auto)
res_immunoClust_all_Mosmann_auto <- data.frame(label = clus_immunoClust_all_Mosmann_auto)


write.table(res_immunoClust_Levine_32_auto, 
            file = "../results_auto/immunoClust/immunoClust_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Levine_13_auto, 
            file = "../results_auto/immunoClust/immunoClust_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Nilsson_auto, 
            file = "../results_auto/immunoClust/immunoClust_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Mosmann_auto, 
            file = "../results_auto/immunoClust/immunoClust_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_immunoClust_all_Levine_32_auto, 
            file = "../results_auto/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Levine_13_auto, 
            file = "../results_auto/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Nilsson_auto, 
            file = "../results_auto/immunoClust_all/immunoClust_all_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Mosmann_auto, 
            file = "../results_auto/immunoClust_all/immunoClust_all_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_immunoClust_auto <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_Levine_32_auto["elapsed"], 
  Levine_2015_marrow_13 = runtime_immunoClust_Levine_13_auto["elapsed"], 
  Nilsson_2013_HSC = runtime_immunoClust_Nilsson_auto["elapsed"], 
  Mosmann_2014_activ = runtime_immunoClust_Mosmann_auto["elapsed"], 
  row.names = "runtime"))

runtime_immunoClust_all_auto <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_all_Levine_32_auto["elapsed"], 
  Levine_2015_marrow_13 = runtime_immunoClust_all_Levine_13_auto["elapsed"], 
  Nilsson_2013_HSC = runtime_immunoClust_all_Nilsson_auto["elapsed"], 
  Mosmann_2014_activ = runtime_immunoClust_all_Mosmann_auto["elapsed"], 
  row.names = "runtime"))

write.table(runtime_immunoClust, 
            file = "../results_auto/runtime/runtime_immunoClust.txt", 
            quote = FALSE, sep = "\t")

write.table(runtime_immunoClust_all, 
            file = "../results_auto/runtime/runtime_immunoClust_all.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_auto/session_info/session_info_immunoClust_and_immunoClust_all.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_auto/RData_files/results_immunoClust_and_immunoClust_all.RData")




#################################################################################
### Run immunoClust and immunoClust_all: manually selected number of clusters ###
#################################################################################

# run immunoClust with manual selection of number of clusters
# (note: decreasing the bias argument increases the number of clusters)

set.seed(123)
runtime_immunoClust_Levine_32_manual <- system.time({
  out_immunoClust_Levine_32_manual <- 
    immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32, bias = 0.1)
})

set.seed(123)
runtime_immunoClust_Levine_13_manual <- system.time({
  out_immunoClust_Levine_13_manual <- 
    immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13, bias = 0.1)
})

set.seed(123)
runtime_immunoClust_Nilsson_manual <- system.time({
  out_immunoClust_Nilsson_manual <- 
    immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, bias = 0.1)
})

set.seed(123)
runtime_immunoClust_Mosmann_manual <- system.time({
  out_immunoClust_Mosmann_manual <- 
    immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, bias = 0.1)
})


# run immmunoClust_all (additional step to classify all cells)

set.seed(123)
runtime_immunoClust_all_Levine_32_manual <- system.time({
  out_immunoClust_all_Levine_32_manual <- 
    immunoClust::cell.process(data_Levine_32, parameters = pars_Levine_32, bias = 0.1, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Levine_13_manual <- system.time({
  out_immunoClust_all_Levine_13_manual <- 
    immunoClust::cell.process(data_Levine_13, parameters = pars_Levine_13, bias = 0.1, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Nilsson_manual <- system.time({
  out_immunoClust_all_Nilsson_manual <- 
    immunoClust::cell.process(data_Nilsson, parameters = pars_Nilsson, bias = 0.1, classify.all = TRUE)
})

set.seed(123)
runtime_immunoClust_all_Mosmann_manual <- system.time({
  out_immunoClust_all_Mosmann_manual <- 
    immunoClust::cell.process(data_Mosmann, parameters = pars_Mosmann, bias = 0.1, classify.all = TRUE)
})


# extract cluster labels

# immunoClust

summary(out_immunoClust_Levine_32_manual)  # number of clusters
summary(out_immunoClust_Levine_13_manual)
summary(out_immunoClust_Nilsson_manual)
summary(out_immunoClust_Mosmann_manual)

clus_immunoClust_Levine_32_manual <- out_immunoClust_Levine_32_manual@label
clus_immunoClust_Levine_13_manual <- out_immunoClust_Levine_13_manual@label
clus_immunoClust_Nilsson_manual <- out_immunoClust_Nilsson_manual@label
clus_immunoClust_Mosmann_manual <- out_immunoClust_Mosmann_manual@label


# immunoClust_all

summary(out_immunoClust_all_Levine_32_manual)
summary(out_immunoClust_all_Levine_13_manual)
summary(out_immunoClust_all_Nilsson_manual)
summary(out_immunoClust_all_Mosmann_manual)

clus_immunoClust_all_Levine_32_manual <- out_immunoClust_all_Levine_32_manual@label
clus_immunoClust_all_Levine_13_manual <- out_immunoClust_all_Levine_13_manual@label
clus_immunoClust_all_Nilsson_manual <- out_immunoClust_all_Nilsson_manual@label
clus_immunoClust_all_Mosmann_manual <- out_immunoClust_all_Mosmann_manual@label


# cluster sizes and number of clusters

# immunoClust

table(clus_immunoClust_Levine_32_manual)
table(clus_immunoClust_Levine_13_manual)
table(clus_immunoClust_Nilsson_manual)
table(clus_immunoClust_Mosmann_manual)

length(table(clus_immunoClust_Levine_32_manual))
length(table(clus_immunoClust_Levine_13_manual))
length(table(clus_immunoClust_Nilsson_manual))
length(table(clus_immunoClust_Mosmann_manual))

# immunoClust_all

table(clus_immunoClust_all_Levine_32_manual)
table(clus_immunoClust_all_Levine_13_manual)
table(clus_immunoClust_all_Nilsson_manual)
table(clus_immunoClust_all_Mosmann_manual)

length(table(clus_immunoClust_all_Levine_32_manual))
length(table(clus_immunoClust_all_Levine_13_manual))
length(table(clus_immunoClust_all_Nilsson_manual))
length(table(clus_immunoClust_all_Mosmann_manual))


# plots

# immunoClust

data_transf_Levine_32_manual <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_32_manual, data_Levine_32)
data_transf_Levine_13_manual <- immunoClust::trans.ApplyToData(out_immunoClust_Levine_13_manual, data_Levine_13)
data_transf_Nilsson_manual <- immunoClust::trans.ApplyToData(out_immunoClust_Nilsson_manual, data_Nilsson)
data_transf_Mosmann_manual <- immunoClust::trans.ApplyToData(out_immunoClust_Mosmann_manual, data_Mosmann)

png("../results_manual/immunoClust/plot_immunoClust_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine_32_manual, data_transf_Levine_32_manual, N = 1000)
dev.off()

png("../results_manual/immunoClust/plot_immunoClust_Levine_2015_marrow_13.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Levine_13_manual, data_transf_Levine_13_manual, N = 1000)
dev.off()

png("../results_manual/immunoClust/plot_immunoClust_Nilsson_2013_HSC.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Nilsson_manual, data_transf_Nilsson_manual, N = 1000)
dev.off()

png("../results_manual/immunoClust/plot_immunoClust_Mosmann_2014_activ.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_Mosmann_manual, data_transf_Mosmann_manual, N = 1000)
dev.off()


# immunoClust_all

data_transf_all_Levine_32_manual <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_32_manual, data_Levine_32)
data_transf_all_Levine_13_manual <- immunoClust::trans.ApplyToData(out_immunoClust_all_Levine_13_manual, data_Levine_13)
data_transf_all_Nilsson_manual <- immunoClust::trans.ApplyToData(out_immunoClust_all_Nilsson_manual, data_Nilsson)
data_transf_all_Mosmann_manual <- immunoClust::trans.ApplyToData(out_immunoClust_all_Mosmann_manual, data_Mosmann)

png("../results_manual/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_32.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine_32_manual, data_transf_all_Levine_32_manual, N = 1000)
dev.off()

png("../results_manual/immunoClust_all/plot_immunoClust_all_Levine_2015_marrow_13.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Levine_13_manual, data_transf_all_Levine_13_manual, N = 1000)
dev.off()

png("../results_manual/immunoClust_all/plot_immunoClust_all_Nilsson_2013_HSC.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Nilsson_manual, data_transf_all_Nilsson_manual, N = 1000)
dev.off()

png("../results_manual/immunoClust_all/plot_immunoClust_all_Mosmann_2014_activ.png", width = 1000, height = 1000)
immunoClust::splom(out_immunoClust_all_Mosmann_manual, data_transf_all_Mosmann_manual, N = 1000)
dev.off()


# save cluster labels

res_immunoClust_Levine_32_manual <- data.frame(label = clus_immunoClust_Levine_32_manual)
res_immunoClust_Levine_13_manual <- data.frame(label = clus_immunoClust_Levine_13_manual)
res_immunoClust_Nilsson_manual <- data.frame(label = clus_immunoClust_Nilsson_manual)
res_immunoClust_Mosmann_manual <- data.frame(label = clus_immunoClust_Mosmann_manual)

res_immunoClust_all_Levine_32_manual <- data.frame(label = clus_immunoClust_all_Levine_32_manual)
res_immunoClust_all_Levine_13_manual <- data.frame(label = clus_immunoClust_all_Levine_13_manual)
res_immunoClust_all_Nilsson_manual <- data.frame(label = clus_immunoClust_all_Nilsson_manual)
res_immunoClust_all_Mosmann_manual <- data.frame(label = clus_immunoClust_all_Mosmann_manual)


write.table(res_immunoClust_Levine_32_manual, 
            file = "../results_manual/immunoClust/immunoClust_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Levine_13_manual, 
            file = "../results_manual/immunoClust/immunoClust_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Nilsson_manual, 
            file = "../results_manual/immunoClust/immunoClust_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_Mosmann_manual, 
            file = "../results_manual/immunoClust/immunoClust_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(res_immunoClust_all_Levine_32_manual, 
            file = "../results_manual/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_32.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Levine_13_manual, 
            file = "../results_manual/immunoClust_all/immunoClust_all_labels_Levine_2015_marrow_13.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Nilsson_manual, 
            file = "../results_manual/immunoClust_all/immunoClust_all_labels_Nilsson_2013_HSC.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_immunoClust_all_Mosmann_manual, 
            file = "../results_manual/immunoClust_all/immunoClust_all_labels_Mosmann_2014_activ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# save runtime

runtime_immunoClust_manual <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_Levine_32_manual["elapsed"], 
  Levine_2015_marrow_13 = runtime_immunoClust_Levine_13_manual["elapsed"], 
  Nilsson_2013_HSC = runtime_immunoClust_Nilsson_manual["elapsed"], 
  Mosmann_2014_activ = runtime_immunoClust_Mosmann_manual["elapsed"], 
  row.names = "runtime"))

runtime_immunoClust_all_manual <- t(data.frame(
  Levine_2015_marrow_32 = runtime_immunoClust_all_Levine_32_manual["elapsed"], 
  Levine_2015_marrow_13 = runtime_immunoClust_all_Levine_13_manual["elapsed"], 
  Nilsson_2013_HSC = runtime_immunoClust_all_Nilsson_manual["elapsed"], 
  Mosmann_2014_activ = runtime_immunoClust_all_Mosmann_manual["elapsed"], 
  row.names = "runtime"))

write.table(runtime_immunoClust, 
            file = "../results_manual/runtime/runtime_immunoClust.txt", 
            quote = FALSE, sep = "\t")

write.table(runtime_immunoClust_all, 
            file = "../results_manual/runtime/runtime_immunoClust_all.txt", 
            quote = FALSE, sep = "\t")


# save session information

sink(file = "../results_manual/session_info/session_info_immunoClust_and_immunoClust_all.txt")
sessionInfo()
sink()


# save R objects

save.image(file = "../results_manual/RData_files/results_immunoClust_and_immunoClust_all.RData")


