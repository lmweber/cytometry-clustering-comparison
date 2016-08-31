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

file_Mosmann_rare = DATA_DIR + "/Mosmann_rare/data/Mosmann_rare.txt"

# indices of protein marker columns (note Python indices start at 0)
marker_cols_Mosmann_rare = list(range(6, 9)) + list(range(10, 21))

data_Mosmann_rare = numpy.loadtxt(fname = file_Mosmann_rare, delimiter = '\t', skiprows = 1, usecols = marker_cols_Mosmann_rare)

# data_Mosmann_rare.shape



######################
### Run PhenoGraph ###
######################

# note: uses maximum number of available cores

import phenograph

communities_Mosmann_rare, graph_Mosmann_rare, Q_Mosmann_rare = phenograph.cluster(data_Mosmann_rare)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results_stability_random_starts/PhenoGraph"

file_out_Mosmann_rare = OUT_DIR + "/python_out_Mosmann_rare.txt"

numpy.savetxt(fname = file_out_Mosmann_rare, X = communities_Mosmann_rare, fmt = '%i', delimiter = '\t')


