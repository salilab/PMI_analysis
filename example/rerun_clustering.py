import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################
nproc = 10
top_dir =  sys.argv[1] 
analys_dir = top_dir+'/analys/'

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')

################################
# Re-run clustering with new 
# restraint and/or parameters
################################


# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

AT.read_models_info()

AT.do_hdbscan_clustering(['EV_sum', 'XLs_sum', 'Psi_vals_0.01', 'Psi_vals_0.1'],
                         min_cluster_size=500,
                         min_samples=1,
                         skip=2)

exit()


