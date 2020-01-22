import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# Change this to the location of your PMI_analysis folder
sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

c = -1 # Cluster we wish to extract
nproc = 4
top_dir =  sys.argv[1] 
analys_dir = top_dir+'/analys/'

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')

################################
# Extract frames
################################

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Point to the selected_models file
HA = AT.get_models_to_extract('analys/selected_models_A_cluster'+str(c)+'_random.csv')
HB = AT.get_models_to_extract('analys/selected_models_B_cluster'+str(c)+'_random.csv')

rmf_file_out_A = 'A_models_clust'+str(c)+'.rmf3'
rmf_file_out_B = 'B_models_clust'+str(c)+'.rmf3'

# Arguments for do_extract_models_single_rmf:
# HA :: Dataframe object from AT.get_models_to_extract()
# file_out :: The output rmf file
AT.do_extract_models_single_rmf(HA, 
                            rmf_file_out_A, # RMF file outputted for Sample A
                            top_dir,        # Top directory containing the PMI output folders
                            analys_dir,     # Analysis directory to write RMF and score files
    scores_prefix = "A_models_clust"+str(c))  # Prefix for the scores file

AT.do_extract_models_single_rmf(HB, rmf_file_out_B, top_dir, analys_dir, scores_prefix = "B_models_clust"+str(c))

exit()


