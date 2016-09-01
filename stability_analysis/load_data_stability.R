#########################################################################################
# Stability analysis: R script to load data
#
# Lukas Weber, September 2016
#########################################################################################


library(flowCore)

DATA_DIR <- "../../../benchmark_data_sets"


# files
files <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare.fcs")
)

# files: non-transformed
files_notransform <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim_notransform.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare_notransform.fcs")
)

is_rare <- c(FALSE, TRUE)  ## Levine_32dim and Mosmann_rare only


# load data

data <- data_notransform <- vector("list", length(files))
names(data) <- names(data_notransform) <- names(files)

for (i in 1:length(data)) {
  data[[i]] <- flowCore::exprs(flowCore::read.FCS(files[[i]], transformation = FALSE, truncate_max_range = FALSE))
  data_notransform[[i]] <- flowCore::exprs(flowCore::read.FCS(files_notransform[[i]], transformation = FALSE, truncate_max_range = FALSE))
}

sapply(data, dim)
sapply(data_notransform, dim)


