#########################################################################################
# R script to run FlowSOM_pre_meta and FlowSOM
#
# Lukas Weber, July 2016
#########################################################################################


library(flowCore)
library(FlowSOM)




#################
### LOAD DATA ###
#################

# filenames

DATA_DIR <- "../../benchmark_data_sets"

files <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim/data/Levine_32dim.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim/data/Levine_13dim.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik/data/Samusik_01.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik/data/Samusik_all.fcs"), 
  Nilsson_rare = file.path(DATA_DIR, "Nilsson_rare/data/Nilsson_rare.fcs"), 
  Mosmann_rare = file.path(DATA_DIR, "Mosmann_rare/data/Mosmann_rare.fcs"), 
  FlowCAP_ND   = file.path(DATA_DIR, "FlowCAP_ND/data/FlowCAP_ND.fcs"), 
  FlowCAP_WNV  = file.path(DATA_DIR, "FlowCAP_WNV/data/FlowCAP_WNV.fcs")
)

# FlowCAP data sets are treated separately since they require clustering algorithms to be
# run individually for each sample

is_FlowCAP <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)


# FlowSOM requires input data as flowFrame objects

data <- vector("list", length(files))
names(data) <- names(files)

for (i in 1:length(data)) {
  f <- files[[i]]
  
  if (!is_FlowCAP[i]) {
    data[[i]] <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
    
  } else {
    smp <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
    smp <- smp[, "sample"]
    d <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
    data[[i]] <- split(d, smp)
  }
}

head(data[[1]])
head(data[[8]][[1]])

sapply(data, length)

sapply(data[!is_FlowCAP], dim)

sapply(data[is_FlowCAP], function(d) {
  sapply(d, function(d2) {
    dim(d2)
  })
})


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




##########################################################
### Run FlowSOM_pre_meta: automatic number of clusters ###
##########################################################

# run FlowSOM_pre_meta with default number of clusters for all data sets

# default is 10x10 grid, i.e. 100 clusters

seed <- 1000
out <- runtimes <- vector("list", length(data))
names(out) <- names(runtimes) <- names(data)

for (i in 1:length(data)) {
  
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      fSOM <- FlowSOM::ReadInput(data[[i]], transform = FALSE, scale = FALSE)
      fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols[[i]])
      fSOM <- FlowSOM::BuildMST(fSOM)
    })
    out[[i]] <- fSOM
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        fSOM <- FlowSOM::ReadInput(data[[i]][[j]], transform = FALSE, scale = FALSE)
        fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols[[i]])
        fSOM <- FlowSOM::BuildMST(fSOM)
      })
      out[[i]][[j]] <- fSOM
    }
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# store output and runtimes for meta-clustering step below
out_pre_meta_auto <- out
runtimes_pre_meta_auto <- runtimes

# example of FlowSOM plots (one data set only)
FlowSOM::PlotStars(out[[1]])

# example showing how to extract cluster labels (one data set only)
str(out[[1]]$map)
head(out[[1]]$map$mapping)
dim(out[[1]]$map$mapping)

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]]$map$mapping[, 1]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- lapply(out[[i]], function(o) o$map$mapping[, 1])
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1.3
    names_i <- rep(names(clus_list_i), times = sapply(clus_list_i, length))
    clus_collapse_i <- unlist(clus_list_i, use.names = FALSE)
    clus[[i]] <- paste(names_i, clus_collapse_i, sep = "_")
  }
}

sapply(clus, length)

# cluster sizes and number of clusters
# (for FlowCAP data sets, total no. of clusters = no. samples * no. clusters per sample)
table(clus[[1]])
sapply(clus, function(cl) length(table(cl)))

# save cluster labels
files_labels <- paste0("../results_auto/FlowSOM_pre_meta/FlowSOM_pre_meta_labels_", 
                       names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../results_auto/runtime/runtime_FlowSOM_pre_meta.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../results_auto/session_info/session_info_FlowSOM_pre_meta.txt")
sessionInfo()
sink()

# save R objects
save.image(file = "../results_auto/RData_files/results_FlowSOM_pre_meta.RData")




##################################################################
### Run FlowSOM_pre_meta: manually selected number of clusters ###
##################################################################

# run FlowSOM_pre_meta with manually selected number of clusters

# grid size 20x20 (400 clusters) for Mosmann_rare (data set with very rare population);
# and grid size 10x10 (100 clusters, i.e. default) for all other data sets

# grid sizes
grid_size <- list(
  Levine_32dim = 10, 
  Levine_13dim = 10, 
  Samusik_01   = 10, 
  Samusik_all  = 10, 
  Nilsson_rare = 10, 
  Mosmann_rare = 20, 
  FlowCAP_ND   = 10, 
  FlowCAP_WNV  = 10
)

seed <- 1000
out <- runtimes <- vector("list", length(data))
names(out) <- names(runtimes) <- names(data)

for (i in 1:length(data)) {
  
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      fSOM <- FlowSOM::ReadInput(data[[i]], transform = FALSE, scale = FALSE)
      fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols[[i]], 
                                xdim = grid_size[[i]], ydim = grid_size[[i]])
      fSOM <- FlowSOM::BuildMST(fSOM)
    })
    out[[i]] <- fSOM
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        fSOM <- FlowSOM::ReadInput(data[[i]][[j]], transform = FALSE, scale = FALSE)
        fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols[[i]], 
                                  xdim = grid_size[[i]], ydim = grid_size[[i]])
        fSOM <- FlowSOM::BuildMST(fSOM)
      })
      out[[i]][[j]] <- fSOM
    }
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# store output and runtimes for meta-clustering step below
out_pre_meta_manual <- out
runtimes_pre_meta_manual <- runtimes

# example of FlowSOM plots (one data set only)
FlowSOM::PlotStars(out[[1]])

# example showing how to extract cluster labels (one data set only)
str(out[[1]]$map)
head(out[[1]]$map$mapping)
dim(out[[1]]$map$mapping)

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]]$map$mapping[, 1]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- lapply(out[[i]], function(o) o$map$mapping[, 1])
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1.3
    names_i <- rep(names(clus_list_i), times = sapply(clus_list_i, length))
    clus_collapse_i <- unlist(clus_list_i, use.names = FALSE)
    clus[[i]] <- paste(names_i, clus_collapse_i, sep = "_")
  }
}

sapply(clus, length)

# cluster sizes and number of clusters
# (for FlowCAP data sets, total no. of clusters = no. samples * no. clusters per sample)
table(clus[[1]])
sapply(clus, function(cl) length(table(cl)))

# save cluster labels
files_labels <- paste0("../results_manual/FlowSOM_pre_meta/FlowSOM_pre_meta_labels_", 
                       names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../results_manual/runtime/runtime_FlowSOM_pre_meta.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../results_manual/session_info/session_info_FlowSOM_pre_meta.txt")
sessionInfo()
sink()

# save R objects
save.image(file = "../results_manual/RData_files/results_FlowSOM_pre_meta.RData")




###################################################################################
### Run FlowSOM (additional meta-clustering step): automatic number of clusters ###
###################################################################################

# run FlowSOM (additional meta-clustering step) with automatic selection of number of clusters

# using results from above (stored in object "out_pre_meta_auto")

seed <- 1000
out <- runtimes <- vector("list", length(out_pre_meta_auto))
names(out) <- names(runtimes) <- names(out_pre_meta_auto)

for (i in 1:length(out_pre_meta_auto)) {
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      meta <- FlowSOM::MetaClustering(out_pre_meta_auto[[i]]$map$codes, method = "metaClustering_consensus")
    })
    out[[i]] <- meta
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        meta <- FlowSOM::MetaClustering(out_pre_meta_auto[[i]][[j]]$map$codes, method = "metaClustering_consensus")
      })
      out[[i]][[j]] <- meta
    }
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# combine runtimes
for (i in 1:length(runtimes)) {
  runtimes[[i]] <- runtimes_pre_meta_auto[[i]] + runtimes[[i]]
}

# check cluster labels (one data set only)
out[[1]]
out[[8]][[1]]

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][out_pre_meta_auto[[i]]$map$mapping[, 1]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(out_pre_meta_auto[[i]]))
    names(clus_list_i) <- names(out_pre_meta_auto[[i]])
    for (j in 1:length(clus_list_i)) {
      clus_list_i[[j]] <- out[[i]][[j]][out_pre_meta_auto[[i]][[j]]$map$mapping[, 1]]
    }
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1.3
    names_i <- rep(names(clus_list_i), times = sapply(clus_list_i, length))
    clus_collapse_i <- unlist(clus_list_i, use.names = FALSE)
    clus[[i]] <- paste(names_i, clus_collapse_i, sep = "_")
  }
}

sapply(clus, length)

# cluster sizes and number of clusters
# (for FlowCAP data sets, total no. of clusters = no. samples * no. clusters per sample)
table(clus[[1]])
sapply(clus, function(cl) length(table(cl)))

# save cluster labels
files_labels <- paste0("../results_auto/FlowSOM/FlowSOM_labels_", names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../results_auto/runtime/runtime_FlowSOM.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../results_auto/session_info/session_info_FlowSOM.txt")
sessionInfo()
sink()

# save R objects
save.image(file = "../results_auto/RData_files/results_FlowSOM.RData")




###########################################################################################
### Run FlowSOM (additional meta-clustering step): manually selected number of clusters ###
###########################################################################################

# run FlowSOM (additional meta-clustering step) with manually selected number of clusters

# using results from above (stored in object "out_pre_meta_manual")

# number of clusters k
k <- list(
  Levine_32dim = 40, 
  Levine_13dim = 40, 
  Samusik_01   = 40, 
  Samusik_all  = 40, 
  Nilsson_rare = 40, 
  Mosmann_rare = 40, 
  FlowCAP_ND   = 7, 
  FlowCAP_WNV  = 4
)

seed <- 1000
out <- runtimes <- vector("list", length(out_pre_meta_manual))
names(out) <- names(runtimes) <- names(out_pre_meta_manual)

for (i in 1:length(out_pre_meta_manual)) {
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      meta <- FlowSOM::metaClustering_consensus(out_pre_meta_manual[[i]]$map$codes, k = k[[i]])
    })
    out[[i]] <- meta
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        meta <- FlowSOM::metaClustering_consensus(out_pre_meta_manual[[i]][[j]]$map$codes, k = k[[i]])
      })
      out[[i]][[j]] <- meta
    }
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# combine runtimes
for (i in 1:length(runtimes)) {
  runtimes[[i]] <- runtimes_pre_meta_manual[[i]] + runtimes[[i]]
}

# check cluster labels (one data set only)
out[[1]]
out[[8]][[1]]

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][out_pre_meta_manual[[i]]$map$mapping[, 1]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(out_pre_meta_manual[[i]]))
    names(clus_list_i) <- names(out_pre_meta_manual[[i]])
    for (j in 1:length(clus_list_i)) {
      clus_list_i[[j]] <- out[[i]][[j]][out_pre_meta_manual[[i]][[j]]$map$mapping[, 1]]
    }
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1.3
    names_i <- rep(names(clus_list_i), times = sapply(clus_list_i, length))
    clus_collapse_i <- unlist(clus_list_i, use.names = FALSE)
    clus[[i]] <- paste(names_i, clus_collapse_i, sep = "_")
  }
}

sapply(clus, length)

# cluster sizes and number of clusters
# (for FlowCAP data sets, total no. of clusters = no. samples * no. clusters per sample)
table(clus[[1]])
sapply(clus, function(cl) length(table(cl)))

# save cluster labels
files_labels <- paste0("../results_manual/FlowSOM/FlowSOM_labels_", names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../results_manual/runtime/runtime_FlowSOM.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../results_manual/session_info/session_info_FlowSOM.txt")
sessionInfo()
sink()

# save R objects
save.image(file = "../results_manual/RData_files/results_FlowSOM.RData")



