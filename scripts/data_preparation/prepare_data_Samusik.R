#########################################################################################
# R script to prepare benchmark data sets Samusik_01 to Samusik_10; Samusik_all
#
# This is a 39-dimensional mass cytometry (CyTOF) data set, consisting of 10 replicate 
# bone marrow samples from C57BL/6J mice (samples from 10 different mice). Samusik_01 to 
# Samusik_10 contain data from each sample individually, and Samusik_all contains data
# from all 10 samples. Dimensions are 39 cell-surface markers and transcription factors;
# manually gated cell population labels are available for 24 immune cell populations.
#
# This R script pre-processes the data sets, adds manually gated cell population labels, 
# and exports them in .txt and .fcs formats.
#
# Source: Samusik et al. (2016), "Automated mapping of phenotype space with single-cell 
# data", Nature Methods, 13(6).
#
# Link to paper: http://www.nature.com/nmeth/journal/v13/n6/full/nmeth.3863.html
# Link to data (.zip file): "https://web.stanford.edu/~samusik/Panorama BM 1-10.zip"
#
# Lukas Weber, July 2016
#########################################################################################


# load packages

library(flowCore)  # from Bioconductor




#################
### LOAD DATA ###
#################

# see above for link to download files


# ----------------------
# load population labels
# ----------------------

# load file
file_gating <- "raw_data/Panorama BM 1-10/population_assignments.txt"
data_gating <- read.table(file_gating, header = FALSE, stringsAsFactors = FALSE, sep = "\t")


# remove line with error (line is partially cut off) and reset row names
data_gating[164401, ]
data_gating <- data_gating[-164401, ]
rownames(data_gating) <- NULL

head(data_gating)
tail(data_gating)
dim(data_gating)


# extract sample numbers
sample <- as.factor(sapply(strsplit(data_gating[, 1], split = " "), function(s) s[1]))
str(sample)
length(sample)
levels(sample)


# extract event (cell) numbers
event <- sapply(strsplit(data_gating[, 1], split = " "), function(s) s[3])
event <- as.numeric(event)
head(event)
str(event)

# add 1 to event numbers (event numbers are provided as index-0, but R-based row numbers
# in the FCS files are index-1)
event <- event + 1
head(event)
str(event)

all(event == floor(event))  # check: all integers
sum(event == 1)  # check: multiple events no. 1 (but not necessarily one for every sample 
                 # due to unassigned cells)

# split event numbers into one data frame per sample
event <- split(event, sample)

str(event)


# extract population names
population <- data_gating[, 2]

head(population)
str(population)

# convert population names to factor (required to write FCS file)
population <- as.factor(population)
levels(population)

table(population)

# save file of population names
write.table(
  data.frame(population = 1:length(levels(population)), name = levels(population)), 
  file = "data/population_names.txt", 
  quote = FALSE, 
  row.names = FALSE)

# split population names into one data frame per sample
population <- split(population, sample)

str(population)



# --------------------------------------------------
# load data from FCS files and add population labels
# --------------------------------------------------

# note: NA = unassigned cells (cells not assigned to any population by manual gating)


files <- list.files("raw_data/Panorama BM 1-10", pattern = "\\.fcs$", full.names = TRUE)
files

marker_names <- flowCore::parameters(flowCore::read.FCS(files[1]))$desc

data <- list()

for (i in 1:length(files)) {
  
  data_i <- flowCore::exprs(flowCore::read.FCS(files[i], transformation = FALSE))
  names_i <- gsub("^.*/|\\.fcs$", "", files[i])
  
  colnames(data_i) <- marker_names
  
  # sample number
  data_i <- cbind(data_i, sample = i)
  
  # event (cell) number
  data_i <- cbind(data_i, event = 1:nrow(data_i))
  
  # population labels (NA = unassigned)
  data_i <- cbind(data_i, label = NA)
  data_i[, "label"][event[[i]]] <- population[[i]]
  
  data[[i]] <- data_i
  names(data)[i] <- names_i
}


# number of cells per sample (check against Supp. Table. S2 Excel spreadsheet)
n_cells <- sapply(data, nrow)
n_cells

# number of assigned cells per sample
n_assigned <- sapply(event, length)
n_assigned

# fraction assigned
frac_assigned <- n_assigned / n_cells
frac_assigned

# total number of cells
sum(n_cells)


# number of columns per file (includes some non-protein columns)
sapply(data, ncol)


head(data[[1]])
tail(data[[1]])
summary(data[[1]])




#########################
### ARCSINH TRANSFORM ###
#########################

# arcsinh transform
# using scale factor 5 for CyTOF data (see Bendall et al. 2011, Supp. Fig. S2)

data_notransform <- data

asinh_scale <- 5

cols_to_scale <- 9:47  # from meta-data in FCS files (use read.FCS() to access)
length(cols_to_scale)  # check

for (i in 1:length(data)) {
  data[[i]][, cols_to_scale] <- asinh(data[[i]][, cols_to_scale] / asinh_scale)
}

summary(data[[1]])




###################
### EXPORT DATA ###
###################

# export data for each sample individually, as well as combined file (data_all)


# combine all samples into a single data frame

data_all <- do.call(rbind, data)
data_all_notransform <- do.call(rbind, data_notransform)

dim(data_all)
head(data_all)
tail(data_all)


# filenames

files_txt <- sprintf("data/Samusik_%02d.txt", 1:10)
files_txt_notransform <- sprintf("data/Samusik_%02d_notransform.txt", 1:10)

files_fcs <- sprintf("data/Samusik_%02d.fcs", 1:10)
files_fcs_notransform <- sprintf("data/Samusik_%02d_notransform.fcs", 1:10)


# save data files in .txt format (takes ~5 min)

for (i in 1:length(data)) {
  write.table(data[[i]], file = files_txt[i], quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(data_notransform[[i]], file = files_txt_notransform[i], quote = FALSE, sep = "\t", row.names = FALSE)
}

write.table(data_all, file = "data/Samusik_all.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data_all_notransform, file = "data/Samusik_all_notransform.txt", quote = FALSE, sep = "\t", row.names = FALSE)


# save data files in .fcs format

for (i in 1:length(data)) {
  flowCore::write.FCS(flowCore::flowFrame(data[[i]]), filename = files_fcs[i])
  flowCore::write.FCS(flowCore::flowFrame(data_notransform[[i]]), filename = files_fcs_notransform[i])
}

flowCore::write.FCS(flowCore::flowFrame(data_all), filename = "data/Samusik_all.fcs")
flowCore::write.FCS(flowCore::flowFrame(data_all_notransform), filename = "data/Samusik_all_notransform.fcs")


