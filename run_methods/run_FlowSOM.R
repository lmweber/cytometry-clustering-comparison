#########################################################################################
# R script to run FlowSOM_pre and FlowSOM
#
# Lukas Weber, August 2016
#########################################################################################


library(flowCore)
library(FlowSOM)




#################
### LOAD DATA ###
#################

# filenames

DATA_DIR <- "../../../benchmark_data_sets"

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


# load data files: FlowSOM requires flowFrame objects

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
    data[[i]] <- flowCore::split(d, smp)
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




#####################################################
### Run FlowSOM_pre: automatic number of clusters ###
#####################################################

# run FlowSOM_pre with default number of clusters for all data sets

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
    cat("data set", names(data[i]), ": run complete\n")
    
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
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# store output and runtimes for meta-clustering step below
out_pre_auto <- out
runtimes_pre_auto <- runtimes

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
    names(clus_list_i) <- names(out[[i]])
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1_3
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
files_labels <- paste0("../../results_auto/FlowSOM_pre/FlowSOM_pre_labels_", 
                       names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../../results_auto/runtimes/runtime_FlowSOM_pre.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../../results_auto/session_info/session_info_FlowSOM_pre.txt")
print(sessionInfo())
sink()

cat("FlowSOM_pre automatic : all runs complete\n")




#############################################################
### Run FlowSOM_pre: manually selected number of clusters ###
#############################################################

# run FlowSOM_pre with manually selected number of clusters

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
    cat("data set", names(data[i]), ": run complete\n")
    
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
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# store output and runtimes for meta-clustering step below
out_pre_manual <- out
runtimes_pre_manual <- runtimes

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
    names(clus_list_i) <- names(out[[i]])
    
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
files_labels <- paste0("../../results_manual/FlowSOM_pre/FlowSOM_pre_labels_", 
                       names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../../results_manual/runtimes/runtime_FlowSOM_pre.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../../results_manual/session_info/session_info_FlowSOM_pre.txt")
print(sessionInfo())
sink()

cat("FlowSOM_pre manual : all runs complete\n")




###################################################################################
### Run FlowSOM (additional meta-clustering step): automatic number of clusters ###
###################################################################################

# run FlowSOM (additional meta-clustering step) with automatic selection of number of clusters

# using results from above (stored in object "out_pre_auto")

seed <- 1000
out <- runtimes <- vector("list", length(out_pre_auto))
names(out) <- names(runtimes) <- names(out_pre_auto)

for (i in 1:length(out_pre_auto)) {
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      meta <- FlowSOM::MetaClustering(out_pre_auto[[i]]$map$codes, method = "metaClustering_consensus")
    })
    out[[i]] <- meta
    cat("data set", names(data[i]), ": run complete\n")
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        meta <- FlowSOM::MetaClustering(out_pre_auto[[i]][[j]]$map$codes, method = "metaClustering_consensus")
      })
      out[[i]][[j]] <- meta
    }
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# combine runtimes
for (i in 1:length(runtimes)) {
  runtimes[[i]] <- runtimes_pre_auto[[i]] + runtimes[[i]]
}

# check cluster labels (one data set only)
out[[1]]
out[[8]][[1]]

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][out_pre_auto[[i]]$map$mapping[, 1]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(out_pre_auto[[i]]))
    names(clus_list_i) <- names(out_pre_auto[[i]])
    for (j in 1:length(clus_list_i)) {
      clus_list_i[[j]] <- out[[i]][[j]][out_pre_auto[[i]][[j]]$map$mapping[, 1]]
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
files_labels <- paste0("../../results_auto/FlowSOM/FlowSOM_labels_", names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../../results_auto/runtimes/runtime_FlowSOM.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../../results_auto/session_info/session_info_FlowSOM.txt")
print(sessionInfo())
sink()

cat("FlowSOM automatic : all runs complete\n")




###########################################################################################
### Run FlowSOM (additional meta-clustering step): manually selected number of clusters ###
###########################################################################################

# run FlowSOM (additional meta-clustering step) with manually selected number of clusters

# using results from above (stored in object "out_pre_manual")

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
out <- runtimes <- vector("list", length(out_pre_manual))
names(out) <- names(runtimes) <- names(out_pre_manual)

for (i in 1:length(out_pre_manual)) {
  if (!is_FlowCAP[i]) {
    runtimes[[i]] <- system.time({
      # note: In the current version of FlowSOM, the recommended function 
      # FlowSOM::metaClustering_consensus() does not pass along the seed argument 
      # correctly, so results are not reproducible. To get around this, we use the 
      # dependency function ConsensusClusterPlus::ConsensusClusterPlus() instead. 
      # However, this will be fixed in the next update of FlowSOM (version 1.5); after 
      # the update the following (simpler) line of code can be used instead.
      #meta <- FlowSOM::metaClustering_consensus(out_pre_manual[[i]]$map$codes, k = k[[i]], seed = seed)
      
      meta <- suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(t(out_pre_manual[[i]]$map$codes), maxK = k[[i]], seed = seed)
      )
      meta <- meta[[k[[i]]]]$consensusClass
    })
    out[[i]] <- meta
    cat("data set", names(data[i]), ": run complete\n")
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      runtimes[[i]][[j]] <- system.time({
        # note: In the current version of FlowSOM, the recommended function 
        # FlowSOM::metaClustering_consensus() does not pass along the seed argument 
        # correctly, so results are not reproducible. To get around this, we use the 
        # dependency function ConsensusClusterPlus::ConsensusClusterPlus() instead. 
        # However, this will be fixed in the next update of FlowSOM (version 1.5); after 
        # the update the following (simpler) line of code can be used instead.
        #meta <- FlowSOM::metaClustering_consensus(out_pre_manual[[i]][[j]]$map$codes, k = k[[i]], seed = seed)
        
        meta <- suppressMessages(
          ConsensusClusterPlus::ConsensusClusterPlus(t(out_pre_manual[[i]][[j]]$map$codes), maxK = k[[i]], seed = seed)
        )
        meta <- meta[[k[[i]]]]$consensusClass
      })
      out[[i]][[j]] <- meta
    }
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# combine runtimes
for (i in 1:length(runtimes)) {
  runtimes[[i]] <- runtimes_pre_manual[[i]] + runtimes[[i]]
}

# check cluster labels (one data set only)
out[[1]]
out[[8]][[1]]

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][out_pre_manual[[i]]$map$mapping[, 1]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(out_pre_manual[[i]]))
    names(clus_list_i) <- names(out_pre_manual[[i]])
    for (j in 1:length(clus_list_i)) {
      clus_list_i[[j]] <- out[[i]][[j]][out_pre_manual[[i]][[j]]$map$mapping[, 1]]
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
files_labels <- paste0("../../results_manual/FlowSOM/FlowSOM_labels_", names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../../results_manual/runtimes/runtime_FlowSOM.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../../results_manual/session_info/session_info_FlowSOM.txt")
print(sessionInfo())
sink()

cat("FlowSOM manual : all runs complete\n")



