#########################################################################################
# R script to prepare benchmark data set Levine_13dim
# 
# This is a 13-dimensional mass cytometry (CyTOF) data set, consisting of expression
# levels of 13 surface marker proteins. Cluster labels are available for 24 manually
# gated cell populations. Cells are healthy human bone marrow mononuclear cells (BMMCs),
# from 1 individual.
#
# This R script pre-processes the data set, adds manually gated cell population labels, 
# and exports it in .txt and .fcs formats.
#
# Source: "benchmark data set 1" in the following paper:
# Levine et al. (2015), "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like
# Cells that Correlate with Prognosis", Cell, 162, 184-197.
#
# Link to paper: http://www.sciencedirect.com/science/article/pii/S0092867415006376
# Link to data: https://www.cytobank.org/cytobank/experiments/46259 (download the FCS
# files with Actions -> Export -> Download Files -> All FCS Files)
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
# one FCS file per manually gated cluster
# 13 surface markers (dimensions), 24 manually gated populations


# FCS filenames
# "unassigned" cells are those where cluster labels are unavailable

files <- list.files("raw_data", pattern = "\\.fcs$", full.names = TRUE)

files_assigned <- files[-grep("NotGated", files)]
files_unassigned <- files[grep("NotGated", files)]

files_assigned
files_unassigned


# cell population names

files_assigned %>% 
  gsub("^.*Marrow1_", "", .) %>% 
  gsub("\\.fcs$", "", .) %>% 
  gsub(" ", "_", .) -> 
  pop_names

pop_names

df_pop_names <- data.frame(label = 1:length(pop_names), population = pop_names)
df_pop_names


# column names (protein markers)

read.FCS(files_assigned[1], transformation = FALSE, truncate_max_range = FALSE) %>% 
  exprs %>% 
  colnames %>% 
  unname -> 
  col_names

col_names


# load FCS files and add cluster labels ("assigned" cells only)

data <- matrix(nrow = 0, ncol = length(col_names) + 1)

for (i in 1:length(files_assigned)) {
  data_i <- flowCore::exprs(flowCore::read.FCS(files_assigned[i], 
                                               transformation = FALSE, 
                                               truncate_max_range = FALSE))
  
  # cluster labels
  data_i <- cbind(data_i, label = i)
  
  data <- rbind(data, data_i)
}

head(data)
dim(data)  # 81,747 assigned cells, 13 dimensions (plus one column of cluster labels)
table(data[, "label"])  # 24 manually gated clusters


# load FCS file for unassigned cells

data_unassigned <- flowCore::exprs(flowCore::read.FCS(files_unassigned, 
                                                      transformation = FALSE, 
                                                      truncate_max_range = FALSE))

data_unassigned <- cbind(data_unassigned, label = NA)

head(data_unassigned)
dim(data_unassigned)  # 85,297 unassigned cells




#########################
### ARCSINH TRANSFORM ###
#########################

# arcsinh transform
# using scale factor 5 for CyTOF data (see Bendall et al. 2011, Supp. Fig. S2)

data_notransform <- data
data_notransform_unassigned <- data_unassigned

asinh_scale <- 5

cols_to_scale <- 1:13
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

write.table(df_pop_names, file = "data/population_names_Levine_13dim.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# save data files in TXT format

write.table(data_combined, file = "data/Levine_13dim.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_combined_notransform, file = "data/Levine_13dim_notransform.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# save data files in FCS format

flowCore::write.FCS(flowCore::flowFrame(data_combined), filename = "data/Levine_13dim.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_combined_notransform), filename = "data/Levine_13dim_notransform.fcs")


