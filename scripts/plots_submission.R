#########################################################################################
# R script to convert plots to EPS format for journal submission
#
# Lukas M. Weber, March 2016
#########################################################################################

# for journal submission: convert plots to EPS format, max width 7.5in, all fonts as outlines



# Figure 1 (Levine_32 multi-panel)

fn.pdf <- "../plots/Levine_2015_marrow_32/plots_multi_panel_Levine2015marrow32.pdf"
fn.eps <- "../plots_submission/Fig1_Levine32.eps"
fn.tmp.eps <- "../plots_submission/Fig1_Levine32.tmp.eps"

system(paste("pdftops -eps -paperw 2100", fn.pdf, fn.tmp.eps))
system(paste("gs -o", fn.eps, "-sDEVICE=epswrite -dNoOutputFonts -dEPSCrop", fn.tmp.eps))
unlink(fn.tmp.eps)



# Figure 2 (Levine_13 multi-panel)

fn.pdf <- "../plots/Levine_2015_marrow_13/plots_multi_panel_Levine2015marrow13.pdf"
fn.eps <- "../plots_submission/Fig2_Levine13.eps"
fn.tmp.eps <- "../plots_submission/Fig2_Levine13.tmp.eps"

system(paste("pdftops -eps -paperw 2100", fn.pdf, fn.tmp.eps))
system(paste("gs -o", fn.eps, "-sDEVICE=epswrite -dNoOutputFonts -dEPSCrop", fn.tmp.eps))
unlink(fn.tmp.eps)



# Figure 3 (Cluster medians for FlowSOM_meta)

fn.pdf <- "../plots/Levine_2015_marrow_32/cluster_medians/cluster_medians_FlowSOM_meta_Levine2015marrow32.pdf"
fn.eps <- "../plots_submission/Fig3_cluster_medians.eps"
fn.tmp.eps <- "../plots_submission/Fig3_cluster_medians.tmp.eps"

system(paste("pdftops -eps -paperw 1400", fn.pdf, fn.tmp.eps))
system(paste("gs -o", fn.eps, "-sDEVICE=epswrite -dNoOutputFonts -dEPSCrop", fn.tmp.eps))
unlink(fn.tmp.eps)



# Figure 4 (Nilsson multi-panel)

fn.pdf <- "../plots/Nilsson_2013_HSC/plots_multi_panel_Nilsson2013HSC.pdf"
fn.eps <- "../plots_submission/Fig4_Nilsson.eps"
fn.tmp.eps <- "../plots_submission/Fig4_Nilsson.tmp.eps"

system(paste("pdftops -eps -paperw 970", fn.pdf, fn.tmp.eps))
system(paste("gs -o", fn.eps, "-sDEVICE=epswrite -dNoOutputFonts -dEPSCrop", fn.tmp.eps))
unlink(fn.tmp.eps)



# Figure 5 (Mosmann multi-panel)

fn.pdf <- "../plots/Mosmann_2014_activ/plots_multi_panel_Mosmann2014activ.pdf"
fn.eps <- "../plots_submission/Fig5_Mosmann.eps"
fn.tmp.eps <- "../plots_submission/Fig5_Mosmann.tmp.eps"

system(paste("pdftops -eps -paperw 970", fn.pdf, fn.tmp.eps))
system(paste("gs -o", fn.eps, "-sDEVICE=epswrite -dNoOutputFonts -dEPSCrop", fn.tmp.eps))
unlink(fn.tmp.eps)



# Figure 6 (Stability analysis for Nilsson and Mosmann)

fn.pdf <- "../plots/Mosmann_2014_activ/stability_analysis/stability_multi_panel_Nilsson_Mosmann.pdf"
fn.eps <- "../plots_submission/Fig6_stability.eps"
fn.tmp.eps <- "../plots_submission/Fig6_stability.tmp.eps"

system(paste("pdftops -eps -paperw 2100", fn.pdf, fn.tmp.eps))
system(paste("gs -o", fn.eps, "-sDEVICE=epswrite -dNoOutputFonts -dEPSCrop", fn.tmp.eps))
unlink(fn.tmp.eps)



