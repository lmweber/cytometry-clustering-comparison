# How to run X-shift

# Not available in R, so follow steps below instead.
# Instructions also available at: http://web.stanford.edu/~samusik/vortex/

# Requires: Standalone application "VorteX"

# 1. X-shift requires a large amount of memory, so it is best to close all other
# applications and reset your computer before starting.
# 2. Run VorteX application.
# 3. Click on the green "plus" symbol next to Datasets to import FCS file.
# 4. Select markers to import by clicking on them; then click "Next".
# 5. Select parameters; click "Next" and "Finalize Import".
# 6. Click on the green "plus" symbol next to Clustering.
# 7. Select parameters; click "Go!".
# 8. After clustering is complete, select all clustering results in the "Clustering"
# window (bottom-left), right-click, and select Validation -> find Elbow Point for
# Cluster Number. This is parameter K corresponding to the optimal automatic number of
# clusters (i.e. not the number of clusters itself). Alternatively, select solution with
# desired number of clusters shown in column "Num. Clusters".
# 9. Right-click solution with selected number of clusters in the "Clustering" window 
# (bottom-left) and select "Export As CSV". Select "Include profiles" if you need to keep
# the original data columns. Cluster labels are in the "ClusterID" column, and cell IDs 
# in the "ProfileID" column. Note that cell IDs start at 0 instead of 1, so you need to 
# add 1 to analyze them in R.




###################
### SUBSAMPLING ###
###################

# The current version of X-shift does not appear to export the column of true population
# labels correctly after subsampling. Use code below to subsample data sets and save true
# population labels instead.

# Note: use non-transformed data sets, since it is better to let X-shift apply the 
# arcsinh transform internally (otherwise noise threshold is applied to the 
# pre-transformed data, which ends up removing more points than intended).


library(flowCore)

# load data

DATA_DIR <- "../../../benchmark_data_sets"

files <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim_notransform.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim_notransform.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01_notransform.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik/data/Samusik_all_notransform.fcs"), 
  Nilsson_rare = file.path(DATA_DIR, "Nilsson_rare/data/Nilsson_rare_notransform.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare_notransform.fcs")
)

data <- vector("list", length(files))
names(data) <- names(files)

for (i in 1:length(data)) {
  f <- files[[i]]
  data[[i]] <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
}

# subsampling for data sets with excessive runtime (> 6 hrs on laptop)

ix_subsample <- 4
n_sub <- 300000

for (i in ix_subsample) {
  set.seed(123)
  data[[i]] <- data[[i]][sample(1:nrow(data[[i]]), n_sub), ]
  
  # save subsampled data sets in FCS format with population labels
  files_sub_i <- paste0(c("../../results/auto/Xshift/", 
                          "../../results/manual/Xshift/"), 
                        names(data)[i], "_notransform_subsampled.fcs")
  for (f in files_sub_i) {
    flowCore::write.FCS(flowCore::flowFrame(data[[i]]), filename = f)
  }
}



