#########################################################################################
# R script to generate plot of results for FlowSOM using range of values for number of
# clusters k
#
# Lukas Weber, September 2016
#########################################################################################


library(ggplot2)
library(reshape2)

# color-blind friendly palettes (http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cb_pal_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_pal_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




####################
### LOAD RESULTS ###
####################

RES_DIR <- "../../results_range_k"

datasets <- c("Levine_32dim", "Levine_13dim", "Samusik_01", "Samusik_all", "Nilsson_rare", "Mosmann_rare")

res <- vector("list", length(datasets))
names(res) <- datasets

for (i in 1:length(res)) {
  file <- paste0(RES_DIR, "/results_FlowSOM_range_k_", datasets[i], ".txt")
  res[[i]] <- read.table(file, header = TRUE, sep = "\t", comment.char = "")
}




#################
### LINE PLOT ###
#################

# plot mean F1 score or F1 score vs. number of clusters k


# tidy data format for ggplot

res_F1 <- data.frame(k = res[[1]][, "k"], 
                     Levine_32dim = res[["Levine_32dim"]][, "mean_F1"], 
                     Levine_13dim = res[["Levine_13dim"]][, "mean_F1"], 
                     Samusik_01 = res[["Samusik_01"]][, "mean_F1"], 
                     Samusik_all = res[["Samusik_all"]][, "mean_F1"], 
                     Nilsson_rare = res[["Nilsson_rare"]][, "F1"], 
                     Mosmann_rare = res[["Mosmann_rare"]][, "F1"])

res_F1_tidy <- melt(res_F1, id.vars = "k", measure_vars = datasets)


# line plot

pl <- 
  ggplot(res_F1_tidy, 
         aes(x = k, y = value, group = variable)) + 
  geom_point(aes(color = variable)) + 
  geom_line(aes(color = variable)) + 
  geom_vline(xintercept = 40, color = "red", linetype = 1) + 
  scale_color_manual(values = cb_pal_black[c(1:6)]) + 
  ylim(0, 0.82) + 
  scale_x_continuous(breaks = seq(5, 80, by = 5)) + 
  ylab("mean F1 score or F1 score") + 
  xlab("number of clusters k") + 
  ggtitle("FlowSOM") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12), 
        legend.position = "right", 
        legend.key.size = unit(5, "mm"), 
        #legend.position = c(0.80, 0.2), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_rect("white"))

print(pl)

ggplot2::ggsave("../../plots/plots_range_k/results_range_k_FlowSOM_F1.pdf", width = 7.5, height = 5)



