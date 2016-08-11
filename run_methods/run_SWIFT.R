# How to run SWIFT

# Not available in R, so follow steps below instead:

# Requires: Matlab, Statistics Toolbox, Parallel Computing Toolbox

# 1. Copy FCS file into a temporary directory, since results files will also be saved 
# here.
# 2. Run SWIFT graphical interface from Matlab by typing "swift_main" in command window.
# 3. After graphical interface opens, select FCS file to import data. Note that SWIFT
# will automatically perform an arcsinh transform, so the FCS file should not be
# transformed already.
# 4. Enter parameters and click to continue.
# 5. After SWIFT completes, cluster labels will be saved in the file
# "<original_filename>.Cluster_Output.txt" in the input directory. Cluster labels are in 
# the "MergeCluster" column.
