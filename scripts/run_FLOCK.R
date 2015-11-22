# How to run FLOCK

# Not available in R, so follow steps below instead:

# Requires: C source code (note that FLOCK can also be run through the ImmPort online
# analysis platform)

# 1. Download C source code from:
# http://sourceforge.net/projects/immportflock/files/FLOCK_flowCAP-I_code/. Be careful to
# download version 2.0, not 1.0. Compilation instructions are in the README file.
# 2. Copy data file into the same folder. The data file should be in .txt format, and
# contain only columns of transformed protein expression values.
# 3. Run from command line with: "./flock2 ./filename.txt"
# 4. Results (population IDs for each cell) are saved in the file flock_results.txt.
