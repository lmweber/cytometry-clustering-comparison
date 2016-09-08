# How to run SWIFT

# Not available in R, so follow steps below instead:

# Requires: Matlab, Statistics Toolbox, Parallel Computing Toolbox

# 1. Copy FCS file into a temporary directory; results files will be saved in the same
# directory.
# 2. Run SWIFT graphical interface from Matlab by typing "swift_main" in command window.
# 3. After graphical interface opens, select FCS file to import data. Note that SWIFT
# will automatically perform an arcsinh transform, so the FCS file should not be
# transformed already.
# 4. Enter parameters and click to continue.
# 5. After SWIFT completes, cluster labels will be saved in the file
# "<original_filename>.Cluster_Output.txt" in the input directory. Cluster labels are in 
# the "MergeCluster" column.




###################
### SUBSAMPLING ###
###################

# SWIFT requires subsampling for some data sets due to runtime. Use code below to
# subsample and save true population labels.

# note: use non-transformed data files

library(flowCore)

# load data

DATA_DIR <- "../../../benchmark_data_sets"

files <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim_notransform.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim_notransform.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01_notransform.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik/data/Samusik_all_notransform.fcs"), 
  Nilsson_rare = file.path(DATA_DIR, "Nilsson_rare/data/Nilsson_rare_notransform.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare_notransform.fcs"), 
  FlowCAP_ND   = file.path(DATA_DIR, "FlowCAP_ND/data/FlowCAP_ND.fcs"), 
  FlowCAP_WNV  = file.path(DATA_DIR, "FlowCAP_WNV/data/FlowCAP_WNV.fcs")
)

is_FlowCAP <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)

data <- vector("list", length(files))
names(data) <- names(files)

for (i in 1:length(data)) {
  f <- files[[i]]
  
  if (!is_FlowCAP[i]) {
    data[[i]] <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
    
  } else {
    smp <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
    smp <- smp[, "sample"]
    d <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
    d <- flowCore::split(d, smp)
    data[[i]] <- lapply(d, function(s) flowCore::exprs(s))
  }
}

# subsampling for data sets with excessive runtime (> 6 hrs on laptop)

ix_subsample <- c(1, 4)
n_sub <- c(100000, NA, NA, 100000, NA, NA)

for (i in ix_subsample) {
  if (!is_FlowCAP[i]) {
    set.seed(123)
    data[[i]] <- data[[i]][sample(1:nrow(data[[i]]), n_sub[i]), ]
    
    # save subsampled data sets in FCS format with population labels
    files_sub_i <- paste0("../../results/auto/SWIFT/", names(data)[i], "_notransform_subsampled.fcs")
    flowCore::write.FCS(flowCore::flowFrame(data[[i]]), filename = files_sub_i)
  }
}




##############################################
### FIX ERROR FOR FLOW CYTOMETRY DATA SETS ###
##############################################

# Some of the non-marker columns (forward scatter etc) included in the flow cytometry 
# data sets (Nilsson_rare, Mosmann_rare) appear to cause an error while trying to load 
# the data into SWIFT. To fix this, save data sets containing marker columns only.

# indices of protein marker columns

marker_cols <- list(
  Levine_32dim = 5:36, 
  Levine_13dim = 1:13, 
  Samusik_01   = 9:47, 
  Samusik_all  = 9:47, 
  Nilsson_rare = c(5:7, 9:18), 
  Mosmann_rare = c(7:9, 11:21), 
  FlowCAP_ND   = 3:12, 
  FlowCAP_WNV  = 3:8
)
sapply(marker_cols, length)

# subset and save data files: protein marker columns only

for (i in c(5, 6)) {
  data[[i]] <- data[[i]][, marker_cols[[i]]]
  files_i <- paste0("../../results/auto/SWIFT/", names(data)[i], "_notransform_markers_only.fcs")
  flowCore::write.FCS(flowCore::flowFrame(data[[i]]), filename = files_i)
}



