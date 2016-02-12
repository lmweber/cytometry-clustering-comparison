#########################################################################################
# Python3 script to run PhenoGraph for robustness analysis
# 
# This script runs one iteration of PhenoGraph (Python implementation) and saves results
# as text files.
# 
# Lukas M. Weber, February 2016
#########################################################################################



#################
### LOAD DATA ###
#################

import numpy

DATA_DIR = "../../../benchmark_data_sets"

file_Levine_32 = DATA_DIR + "/Levine_2015_marrow_32/data/Levine_2015_marrow_32.txt"

data_Levine_32 = numpy.loadtxt(fname = file_Levine_32, delimiter = '\t', skiprows = 1)


# indices of protein marker columns
# note: Python indices start at 0

marker_cols_Levine_32 = list(range(4, 36))


# subset data

data_Levine_32 = data_Levine_32[:, marker_cols_Levine_32]

# data_Levine_32.shape




######################
### Run PhenoGraph ###
######################

# set n_jobs = 1 to use 1 core only, for comparability with main results

import phenograph

communities_Levine_32, graph_Levine_32, Q_Levine_32 = phenograph.cluster(data_Levine_32, n_jobs = 1)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results/robustness_analysis/PhenoGraph"

file_out_Levine_32 = OUT_DIR + "/python_out_Levine_32.txt"

numpy.savetxt(fname = file_out_Levine_32, X = communities_Levine_32, fmt = '%i', delimiter = '\t')



