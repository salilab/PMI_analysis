import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# Change to the location of your PMI_analysis folder
sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from analysis_trajectories import *

if __name__ == "__main__":
    nproc = 1
    top_dir =  sys.argv[1]
    cluster = sys.argv[2]
    if len(sys.argv) == 4:
        state = sys.argv[3]
    else:
        state = 0

    analysis_dir = os.path.join(top_dir,'analysis/')
    print(analysis_dir)

    # How are the trajectories dir names
    dir_head = 'run_'
    out_dirs = glob.glob(os.path.join(top_dir,f"{dir_head}*/output/"))

    # Load module
    AT = AnalysisTrajectories(out_dirs,
                              dir_name=dir_head,
                              analysis_dir = analysis_dir,
                              nproc=nproc)

    # Point to the selected_models file
    HA = AT.get_models_to_extract(
        f'{analysis_dir}/selected_models_A_cluster{cluster}_detailed.csv')
    HB = AT.get_models_to_extract(
        f'{analysis_dir}/selected_models_B_cluster{cluster}_detailed.csv')

    rmf_file_out_A = f'A_models_cluster{cluster}_{state}.rmf3'
    rmf_file_out_B = f'B_models_cluster{cluster}_{state}.rmf3'

    AT.do_extract_models_single_rmf(HA, 
                                    rmf_file_out_A,
                                    top_dir,       
                                    analysis_dir,     
                                    scores_prefix = f"A_models_cluster{cluster}_{state}",
                                    sel_state = state)  
                         
    AT.do_extract_models_single_rmf(HB,
                                    rmf_file_out_B,
                                    top_dir,
                                    analysis_dir,
                                    scores_prefix = f"B_models_cluster{cluster}_{state}",
                                    sel_state = state)




