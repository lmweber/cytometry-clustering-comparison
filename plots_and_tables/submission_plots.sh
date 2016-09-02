#########################################################################################
# Shell script to copy and rename plots for journal submission
#
# Lukas Weber, September 2016
#########################################################################################


# for journal submission: keep figures in PDF format, but rename to "Fig1.pdf" etc


# Figure 1: Levine_32dim main results
cp ../../plots/Levine_32dim/main_plots/plots_multi_panel_Levine_32dim.pdf ../../plots_submission/Fig1_Levine_32dim.pdf

# Figure 2: Cluster medians (Levine_32dim, FlowSOM)
cp ../../plots/Levine_2015_marrow_32/cluster_medians/cluster_medians_heatmap_FlowSOM_meta_Levine2015marrow32.pdf ../../plots_submission/Fig2_cluster_medians.pdf

# Figure 3: Mosmann_rare main results
cp ../../plots/Mosmann_rare/main_plots/plots_multi_panel_Mosmann_rare.pdf ../../plots_submission/Fig3_Mosmann_rare.pdf

