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

# Create dir
gsms_A_dir = analys_dir+'GSMs_cl0/sample_A'
gsms_B_dir = analys_dir+'GSMs_cl0/sample_B'

AT.create_gsms_dir(gsms_A_dir)
AT.create_gsms_dir(gsms_B_dir)

HA = AT.get_models_to_extract('analys/selected_models_A_cluster0_random.csv')
HB = AT.get_models_to_extract('analys/selected_models_B_cluster0_random.csv')
AT.do_extract_models(HA, 'h1', gsms_A_dir)
AT.do_extract_models(HB, 'h2', gsms_B_dir)

exit()


