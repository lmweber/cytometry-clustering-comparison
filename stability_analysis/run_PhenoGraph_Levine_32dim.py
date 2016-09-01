#########################################################################################
# Python3 script to run PhenoGraph for stability analysis
# 
# Lukas Weber, August 2016
#########################################################################################



#################
### LOAD DATA ###
#################

import numpy

DATA_DIR = "../../../benchmark_data_sets"

file_Levine_32dim = DATA_DIR + "/Levine_32dim/data/Levine_32dim.txt"

# indices of protein marker columns (note Python indices start at 0)
marker_cols_Levine_32dim = list(range(4, 36))

data_Levine_32dim = numpy.loadtxt(fname = file_Levine_32dim, delimiter = '\t', skiprows = 1, usecols = marker_cols_Levine_32dim)

# data_Levine_32dim.shape



######################
### Run PhenoGraph ###
######################

# note: uses maximum number of available cores

import phenograph

communities_Levine_32dim, graph_Levine_32dim, Q_Levine_32dim = phenograph.cluster(data_Levine_32dim)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results_stability_random_starts/PhenoGraph"

file_out_Levine_32dim = OUT_DIR + "/python_out_Levine_32dim.txt"

numpy.savetxt(fname = file_out_Levine_32dim, X = communities_Levine_32dim, fmt = '%i', delimiter = '\t')


