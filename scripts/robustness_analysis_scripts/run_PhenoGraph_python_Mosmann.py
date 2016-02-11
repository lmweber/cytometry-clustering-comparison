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

file_Mosmann = DATA_DIR + "/Mosmann_2014_activ/data/Mosmann_2014_activ.txt"

data_Mosmann = numpy.loadtxt(fname = file_Mosmann, delimiter = '\t', skiprows = 1)


# indices of protein marker columns
# note: Python indices start at 0

marker_cols_Mosmann = list(range(6, 21))


# subset data

data_Mosmann = data_Mosmann[:, marker_cols_Mosmann]

# data_Mosmann.shape




######################
### Run PhenoGraph ###
######################

# set n_jobs = 1 to use 1 core only, for comparability with main results

import phenograph

communities_Mosmann, graph_Mosmann, Q_Mosmann = phenograph.cluster(data_Mosmann, n_jobs = 1)



####################
### SAVE RESULTS ###
####################

# export results as tab-delimited text file

OUT_DIR = "../../results/PhenoGraph"

file_out_Mosmann = OUT_DIR + "/python_out_Mosmann.txt"

numpy.savetxt(fname = file_out_Mosmann, X = communities_Mosmann, fmt = '%i', delimiter = '\t')



