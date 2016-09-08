#########################################################################################
# R script to run flowPeaks
#
# Lukas Weber, July 2016
#########################################################################################


library(flowCore)
library(flowPeaks)




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


# load data files

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


# subset data: protein marker columns only

for (i in 1:length(data)) {
  if (!is_FlowCAP[i]) {
    data[[i]] <- data[[i]][, marker_cols[[i]]]
  } else {
    for (j in 1:length(data[[i]])) {
      data[[i]][[j]] <- data[[i]][[j]][, marker_cols[[i]]]
    }
  }
}

sapply(data[!is_FlowCAP], dim)
sapply(data[is_FlowCAP], function(d) {
  sapply(d, function(d2) {
    dim(d2)
  })
})




###################################################
### Run flowPeaks: automatic number of clusters ###
###################################################

# run flowPeaks with automatic selection of number of clusters

seed <- 123
out <- runtimes <- vector("list", length(data))
names(out) <- names(runtimes) <- names(data)

for (i in 1:length(data)) {
  
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      out[[i]] <- flowPeaks(data[[i]])
    })
    cat("data set", names(data[i]), ": run complete\n")
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        out[[i]][[j]] <- flowPeaks(data[[i]][[j]])
      })
    }
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][["peaks.cluster"]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(data[[i]]))
    names(clus_list_i) <- names(data[[i]])
    for (j in 1:length(data[[i]])) {
      if (!is.null(out[[i]][[j]])) {
        clus_list_i[[j]] <- out[[i]][[j]][["peaks.cluster"]]
      }
    }
    
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
files_labels <- paste0("../../results/auto/flowPeaks/flowPeaks_labels_", 
                       names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../../results/auto/runtimes/runtime_flowPeaks.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../../results/auto/session_info/session_info_flowPeaks.txt")
print(sessionInfo())
sink()

cat("flowPeaks automatic : all runs complete\n")




################################################
### Run flowPeaks: manual number of clusters ###
################################################

# run flowPeaks with manual selection of number of clusters
# (note: smaller "tol" argument increases number of clusters)

seed <- 123
out <- runtimes <- vector("list", length(data))
names(out) <- names(runtimes) <- names(data)

for (i in 1:length(data)) {
  
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      out[[i]] <- flowPeaks(data[[i]], tol = 0.0001)
    })
    cat("data set", names(data[i]), ": run complete\n")
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        out[[i]][[j]] <- flowPeaks(data[[i]][[j]], tol = 0.0001)
      })
    }
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][["peaks.cluster"]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(data[[i]]))
    names(clus_list_i) <- names(data[[i]])
    for (j in 1:length(data[[i]])) {
      if (!is.null(out[[i]][[j]])) {
        clus_list_i[[j]] <- out[[i]][[j]][["peaks.cluster"]]
      }
    }
    
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
files_labels <- paste0("../../results/manual/flowPeaks/flowPeaks_labels_", 
                       names(clus), ".txt")

for (i in 1:length(files_labels)) {
  res_i <- data.frame(label = clus[[i]])
  write.table(res_i, file = files_labels[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save runtimes
runtimes <- lapply(runtimes, function(r) r["elapsed"])
runtimes <- t(as.data.frame(runtimes, row.names = "runtime"))

write.table(runtimes, file = "../../results/manual/runtimes/runtime_flowPeaks.txt", 
            quote = FALSE, sep = "\t")

# save session information
sink(file = "../../results/manual/session_info/session_info_flowPeaks.txt")
print(sessionInfo())
sink()

cat("flowPeaks manual : all runs complete\n")



