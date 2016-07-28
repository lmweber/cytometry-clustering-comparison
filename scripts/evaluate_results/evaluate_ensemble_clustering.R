#########################################################################################
# R script to load and calculate results for ensemble clustering
#
# Lukas M. Weber, March 2016
#########################################################################################


# helper functions
source("helper_match_clusters_and_evaluate.R")
source("helper_match_one_rare_cluster_and_evaluate.R")

# true (manually gated) cluster labels
source("load_results_truth.R")




####################
### LOAD RESULTS ###
####################

# load ensemble clustering results

file_ensemble_Levine_32 <- "../results_ensemble/ensemble_clustering/ensemble_labels_Levine_2015_marrow_32.txt"
file_ensemble_Levine_13 <- "../results_ensemble/ensemble_clustering/ensemble_labels_Levine_2015_marrow_13.txt"
file_ensemble_Nilsson <- "../results_ensemble/ensemble_clustering/ensemble_labels_Nilsson_2013_HSC.txt"
file_ensemble_Mosmann <- "../results_ensemble/ensemble_clustering/ensemble_labels_Mosmann_2014_activ.txt"

data_ensemble_Levine_32 <- read.table(file_ensemble_Levine_32, header = TRUE, sep = "\t")
data_ensemble_Levine_13 <- read.table(file_ensemble_Levine_13, header = TRUE, sep = "\t")
data_ensemble_Nilsson <- read.table(file_ensemble_Nilsson, header = TRUE, sep = "\t")
data_ensemble_Mosmann <- read.table(file_ensemble_Mosmann, header = TRUE, sep = "\t")

cons_clus_Levine_32 <- data_ensemble_Levine_32[, 1]
cons_clus_Levine_13 <- data_ensemble_Levine_13[, 1]
cons_clus_Nilsson <- data_ensemble_Nilsson[, 1]
cons_clus_Mosmann <- data_ensemble_Mosmann[, 1]


# cluster frequency tables

table(cons_clus_Levine_32)
table(cons_clus_Levine_13)
table(cons_clus_Nilsson)
table(cons_clus_Mosmann)


# number of consensus clusters

length(table(cons_clus_Levine_32))
length(table(cons_clus_Levine_13))
length(table(cons_clus_Nilsson))
length(table(cons_clus_Mosmann))


# confusion matrices

table(ensemble = cons_clus_Levine_32, true = clus_truth_Levine_32)
table(ensemble = cons_clus_Levine_13, true = clus_truth_Levine_13)
table(ensemble = cons_clus_Nilsson, true = clus_truth_Nilsson)
table(ensemble = cons_clus_Mosmann, true = clus_truth_Mosmann)




#############################################
### CALCULATE F1 SCORE, PRECISION, RECALL ###
#############################################

# for data sets with multiple populations of interest (Levine_2015_marrow_32, Levine_2015_marrow_13)

# mean F1 score, mean precision, mean recall

res_ensemble_Levine_32 <- helper_match_clusters_and_evaluate(cons_clus_Levine_32, clus_truth_Levine_32)
res_ensemble_Levine_13 <- helper_match_clusters_and_evaluate(cons_clus_Levine_13, clus_truth_Levine_13)

mean(res_ensemble_Levine_32$F1)
mean(res_ensemble_Levine_13$F1)


# for data sets with a single rare cell population of interest (Nilsson_2013_HSC, Mosmann_2014_activ)

# F1 score, precision, recall

res_ensemble_Nilsson <- helper_match_one_rare_cluster_and_evaluate(cons_clus_Nilsson, clus_truth_Nilsson)
res_ensemble_Mosmann <- helper_match_one_rare_cluster_and_evaluate(cons_clus_Mosmann, clus_truth_Mosmann)

res_ensemble_Nilsson$F1
res_ensemble_Mosmann$F1



