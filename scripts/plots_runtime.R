#########################################################################################
# R script to generate plot of runtime
#
# Lukas M. Weber, November 2015
#########################################################################################


library(ggplot2)
library(reshape2)



##########################
### PREPARE DATA FRAME ###
##########################

runtime_Levine <- list()
runtime_Mosmann <- list()


# ---------------------------------------------------------
# runtime for methods that were not available as R packages
# ---------------------------------------------------------

# these methods were timed manually with a stopwatch; results saved in spreadsheet
# parameters_and_runtime.xlsx

runtime_Levine[["ACCENSE"]] <- (5 * 60) + 12
runtime_Mosmann[["ACCENSE"]] <- (4 * 60) + 17

runtime_Levine[["FLOCK"]] <- 30
runtime_Mosmann[["FLOCK"]] <- 65

runtime_Levine[["PhenoGraph"]] <- (6 * 60) + 41
runtime_Mosmann[["PhenoGraph"]] <- (49 * 60) + 32

runtime_Levine[["SWIFT"]] <- 3600 + (55 * 60) + 44
runtime_Mosmann[["SWIFT"]] <- 3600 + (16 * 60) + 50


# -----------------------------
# runtime for all other methods
# -----------------------------

RUNTIME_DIR <- "../results/runtime"

other_names_Levine <- c("DensVM", "flowMeans", "FlowSOM", "FlowSOM_meta", "immunoClust", 
                        "immunoClust_all", "kmeans", "Rclusterpp", "SamSPECTRAL")

other_names_Mosmann <- c("DensVM", "flowMeans", "FlowSOM", "FlowSOM_meta", "immunoClust", 
                         "immunoClust_all", "kmeans", "SamSPECTRAL")  # no Rclusterpp


for (i in 1:length(other_names_Levine)) {
  
  file_i <- paste0(RUNTIME_DIR, "/runtime_", other_names_Levine[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  
  runtime_i <- data_i["Levine_2015_marrow_32", "runtime"]
  
  runtime_Levine[[other_names_Levine[i]]] <- runtime_i
}


for (i in 1:length(other_names_Mosmann)) {
  
  file_i <- paste0(RUNTIME_DIR, "/runtime_", other_names_Mosmann[i], ".txt")
  data_i <- read.table(file_i, header = TRUE, sep = "\t")
  
  runtime_i <- data_i["Mosmann_2014_rare", "runtime"]
  
  runtime_Mosmann[[other_names_Mosmann[i]]] <- runtime_i
}


# ---------------
# order ascending
# ---------------

runtime_Levine_ord <- unlist(runtime_Levine)
runtime_Mosmann_ord <- unlist(runtime_Mosmann)

runtime_Levine_ord <- runtime_Levine_ord[order(runtime_Levine_ord)]
runtime_Mosmann_ord <- runtime_Mosmann_ord[order(runtime_Mosmann_ord)]



######################
### GENERATE PLOTS ###
######################

# tidy data format

runtime_Levine_tidy <- data.frame(value = runtime_Levine_ord)
runtime_Levine_tidy["method"] <- factor(rownames(runtime_Levine_tidy), 
                                        levels = rownames(runtime_Levine_tidy))

runtime_Mosmann_tidy <- data.frame(value = runtime_Mosmann_ord)
runtime_Mosmann_tidy["method"] <- factor(rownames(runtime_Mosmann_tidy), 
                                        levels = rownames(runtime_Mosmann_tidy))

# single or multiple cores used

runtime_Levine_tidy$cores <- "single"
runtime_Levine_tidy[c("SWIFT", "Rclusterpp"), "cores"] <- "multiple"

runtime_Mosmann_tidy$cores <- "single"
runtime_Mosmann_tidy[c("SWIFT"), "cores"] <- "multiple"

runtime_Levine_tidy
runtime_Mosmann_tidy


# bar plots

ggplot(runtime_Levine_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "purple2") + 
  geom_text(aes(label = round(value, 0), y = value + 500), size = 3.5) + 
  ggtitle("Runtime: Levine_2015_marrow_32") + 
  ylim(0, 21000) + 
  ylab("seconds") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/runtime_Levine2015marrow32.pdf", height = 6, width = 5.5)


ggplot(runtime_Mosmann_tidy, aes(x = method, y = value)) + 
  geom_bar(stat = "identity", fill = "purple4") + 
  geom_text(aes(label = round(value, 0), y = value + 500), size = 3.5) + 
  ggtitle("Runtime: Mosmann_2014_rare") + 
  ylim(0, 21000) + 
  ylab("seconds") + 
  theme_bw() + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plots/runtime_Mosmann2014rare.pdf", height = 6, width = 5.25)


