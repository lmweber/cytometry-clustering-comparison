#########################################################################################
# Shell script to copy and rename plots for journal submission
#
# Lukas Weber, September 2016
#########################################################################################


# for journal submission: convert to EPS format, and rename to "Fig1.eps" etc


# Figure 1: Levine_32dim main results

# rename and move PDF file
cp ../../plots/Levine_32dim/main_plots/plots_multi_panel_Levine_32dim.pdf ../../plots/submission/Fig1_Levine_32dim.pdf
# convert to EPS
pdftops -eps -paperw 2100 ../../plots/submission/Fig1_Levine_32dim.pdf ../../plots/submission/Fig1_Levine_32dim.eps


# Figure 2: Cluster medians (Levine_32dim, FlowSOM)

# rename and move PDF file
cp ../../plots/Levine_32dim/cluster_medians/cluster_medians_heatmap_FlowSOM_Levine32dim.pdf ../../plots/submission/Fig2_cluster_medians.pdf
# convert to EPS
pdftops -eps -paperw 1400 ../../plots/submission/Fig2_cluster_medians.pdf ../../plots/submission/Fig2_cluster_medians.eps


# Figure 3: Mosmann_rare main results

# rename and move PDF file
cp ../../plots/Mosmann_rare/main_plots/plots_multi_panel_Mosmann_rare.pdf ../../plots/submission/Fig3_Mosmann_rare.pdf
# convert to EPS
pdftops -eps -paperw 970 ../../plots/submission/Fig3_Mosmann_rare.pdf ../../plots/submission/Fig3_Mosmann_rare.eps


