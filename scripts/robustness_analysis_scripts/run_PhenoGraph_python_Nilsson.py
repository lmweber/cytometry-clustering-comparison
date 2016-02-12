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

file_Nilsson = DATA_DIR + "/Nilsson_2013_HSC/data/Nilsson_2013_HSC.txt"

data_Nilsson = numpy.loadtxt(fname = file_Nilsson, delimiter = '\t', skiprows = 1)


# indices of protein marker columns
# note: Python indices start at 0

marker_cols_Nilsson = list(range(4, 7)) + list(range(8, 18))


# subset data

data_Nilsson = data_Nilsson[:, marker_cols_Nilsson]

# data_Nilsson.shape




######################
### Run PhenoGraph ###
######################

# set n_jobs = 1 to use 1 core only, for comparability with main results

import phenograph

communities_Nilsson, graph_Nilsson, Q_Nilsson = phenograph.cluster(data_Nilsson, n_jobs = 1)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results/robustness_analysis/PhenoGraph"

file_out_Nilsson = OUT_DIR + "/python_out_Nilsson.txt"

numpy.savetxt(fname = file_out_Nilsson, X = communities_Nilsson, fmt = '%i', delimiter = '\t')



