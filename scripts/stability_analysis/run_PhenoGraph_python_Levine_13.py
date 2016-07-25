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

file_Levine_13 = DATA_DIR + "/Levine_2015_marrow_13/data/Levine_2015_marrow_13.txt"

# indices of protein marker columns (note Python indices start at 0)
marker_cols_Levine_13 = list(range(0, 13))

data_Levine_13 = numpy.loadtxt(fname = file_Levine_13, delimiter = '\t', skiprows = 1, usecols = marker_cols_Levine_13)

# data_Levine_13.shape



######################
### Run PhenoGraph ###
######################

# note: tried setting n_jobs = 1 for comparability with main results, but doesn't appear to work

import phenograph

communities_Levine_13, graph_Levine_13, Q_Levine_13 = phenograph.cluster(data_Levine_13)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results_stability_analysis/PhenoGraph"

file_out_Levine_13 = OUT_DIR + "/python_out_Levine_13.txt"

numpy.savetxt(fname = file_out_Levine_13, X = communities_Levine_13, fmt = '%i', delimiter = '\t')


