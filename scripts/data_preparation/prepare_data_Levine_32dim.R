#########################################################################################
# R script to prepare benchmark data set Levine_32dim
# 
# This is a 32-dimensional mass cytometry (CyTOF) data set, consisting of expression
# levels of 32 surface marker proteins. Cluster labels are available for 14 manually
# gated cell populations. Cells are healthy human bone marrow mononuclear cells (BMMCs),
# from 2 individuals.
#
# This R script pre-processes the data set, adds manually gated cell population labels, 
# and exports it in .txt and .fcs formats.
#
# Source: "benchmark data set 2" in the following paper:
# Levine et al. (2015), "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like
# Cells that Correlate with Prognosis", Cell, 162, 184-197.
#
# Link to paper: http://www.sciencedirect.com/science/article/pii/S0092867415006376
# Link to data: https://www.cytobank.org/cytobank/experiments/46102 (download the ZIP 
# file shown under "Exported Files")
# 
# Lukas Weber, July 2016
#########################################################################################


# load packages

library(flowCore)  # from Bioconductor
library(magrittr)  # from CRAN




#################
### LOAD DATA ###
#################

# see above for link to download files
# one FCS file per manually gated cluster, per individual (H1 and H2)
# 32 surface markers (dimensions), 14 manually gated populations, 2 individuals (H1 and H2)


# FCS filenames
# "unassigned" cells are those where cluster labels are unavailable

files <- list.files("raw_data", pattern = "\\.fcs$", full.names = TRUE)

files_assigned <- files[-grep("NotDebrisSinglets", files)]
files_unassigned <- files[grep("NotDebrisSinglets", files)]

files_assigned
files_unassigned


# cell population names

grep("H1\\.fcs$", files_assigned) %>% 
  files_assigned[.] %>% 
  gsub("(_H1\\.fcs$)", "", .) %>% 
  gsub("^.*_", "", .) %>% 
  gsub(" ", "_", .) -> 
  pop_names

pop_names

df_pop_names <- data.frame(label = 1:length(pop_names), population = pop_names)
df_pop_names


# column names (protein markers and others)

read.FCS(files_assigned[1], transformation = FALSE, truncate_max_range = FALSE) %>% 
  exprs %>% 
  colnames %>% 
  unname %>% 
  gsub("\\(.*", "", .) -> 
  col_names

col_names


# vector of labels for individuals H1 and H2

indiv <- rep(NA, length(files_assigned))
indiv[grep("_H1\\.fcs$", files_assigned)] <- 1
indiv[grep("_H2\\.fcs$", files_assigned)] <- 2
indiv

indiv_unassigned <- c(1, 2)
indiv_unassigned



# load FCS files and add cluster labels ("assigned" cells only)

data <- matrix(nrow = 0, ncol = length(col_names) + 2)

for (i in 1:length(files_assigned)) {
  data_i <- flowCore::exprs(flowCore::read.FCS(files_assigned[i], 
                                               transformation = FALSE, 
                                               truncate_max_range = FALSE))
  colnames(data_i) <- col_names
  
  # cluster labels
  data_i <- cbind(data_i, label = ((i - 1) %% length(pop_names)) + 1)
  
  # labels for each individual
  data_i <- cbind(data_i, individual = indiv[i])
  
  data <- rbind(data, data_i)
}

head(data)
dim(data)  # 104,184 assigned cells, 32 dimensions (plus 9 other columns)
table(data[, "label"])  # 14 manually gated clusters
table(data[, "individual"])  # 2 individuals (72,463 and 31,721 assigned cells each)



# load FCS files for unassigned cells

data_unassigned <- matrix(nrow = 0, ncol = length(col_names) + 2)

for (i in 1:length(files_unassigned)) {
  data_i <- flowCore::exprs(flowCore::read.FCS(files_unassigned[i], 
                                               transformation = FALSE, 
                                               truncate_max_range = FALSE))
  colnames(data_i) <- col_names
  
  # cluster labels (NA since unassigned)
  data_i <- cbind(data_i, label = NA)
  
  # labels for each individual
  data_i <- cbind(data_i, individual = indiv_unassigned[i])
  
  data_unassigned <- rbind(data_unassigned, data_i)
}

head(data_unassigned)
dim(data_unassigned)  # 161,443 unassigned cells
table(data_unassigned[, "individual"])  # 2 individuals (118,888 and 42,555 unassigned cells each)




#########################
### ARCSINH TRANSFORM ###
#########################

# arcsinh transform
# using scale factor 5 for CyTOF data (see Bendall et al. 2011, Supp. Fig. S2)

data_notransform <- data
data_notransform_unassigned <- data_unassigned

asinh_scale <- 5

cols_to_scale <- 3:38
data[, cols_to_scale] <- asinh(data[, cols_to_scale] / asinh_scale)
data_unassigned[, cols_to_scale] <- asinh(data_unassigned[, cols_to_scale] / asinh_scale)

summary(data)
summary(data_unassigned)




###################
### EXPORT DATA ###
###################

# combine data frames for assigned and unassigned cells

data_combined <- rbind(data, data_unassigned)
data_combined_notransform <- rbind(data_notransform, data_notransform_unassigned)

dim(data_combined)
dim(data_combined_notransform)

# export cell population names

write.table(df_pop_names, file = "data/population_names_Levine_32dim.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# save data files in TXT format

write.table(data_combined, file = "data/Levine_32dim.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_combined_notransform, file = "data/Levine_32dim_notransform.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# save data files in FCS format

flowCore::write.FCS(flowCore::flowFrame(data_combined), filename = "data/Levine_32dim.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_combined_notransform), filename = "data/Levine_32dim_notransform.fcs")


