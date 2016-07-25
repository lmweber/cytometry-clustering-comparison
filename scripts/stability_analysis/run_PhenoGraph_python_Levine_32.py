#########################################################################################
# Python3 script to run PhenoGraph for stability analysis
# 
# This script runs one iteration of PhenoGraph (Python implementation) and saves results
# as text files.
# 
# Lukas M. Weber, March 2016
#########################################################################################



#################
### LOAD DATA ###
#################

import numpy

DATA_DIR = "../../../benchmark_data_sets"

file_Levine_32 = DATA_DIR + "/Levine_2015_marrow_32/data/Levine_2015_marrow_32.txt"

# indices of protein marker columns (note Python indices start at 0)
marker_cols_Levine_32 = list(range(4, 36))

data_Levine_32 = numpy.loadtxt(fname = file_Levine_32, delimiter = '\t', skiprows = 1, usecols = marker_cols_Levine_32)

# data_Levine_32.shape



######################
### Run PhenoGraph ###
######################

# note: tried setting n_jobs = 1 for comparability with main results, but doesn't appear to work

import phenograph

communities_Levine_32, graph_Levine_32, Q_Levine_32 = phenograph.cluster(data_Levine_32)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results_stability_analysis/PhenoGraph"

file_out_Levine_32 = OUT_DIR + "/python_out_Levine_32.txt"

numpy.savetxt(fname = file_out_Levine_32, X = communities_Levine_32, fmt = '%i', delimiter = '\t')


