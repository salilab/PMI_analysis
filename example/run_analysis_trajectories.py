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
    nproc = 10
    top_dir =  sys.argv[1] 
    analys_dir = os.path.join(top_dir,"analysis/")
    
    # Check if analysis dir exists
    if not os.path.isdir(analys_dir):
        os.makedirs(analys_dir)
        
    # How are the trajectories dir names
    dir_head = 'run_'
    out_dirs = glob.glob(
        os.path.join(top_dir,f"{dir_head}*/output/"))
        
    # Add all possible XLs
    XLs_cutoffs = {'DSSO':30.0}
    
    # Load module
    AT = AnalysisTrajectories(out_dirs,
                              dir_name=dir_head,
                              analysis_dir = analys_dir,
                              nproc=nproc)

    # Define restraints to analyze
    AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs)
    AT.set_analyze_Connectivity_restraint()
    AT.set_analyze_Excluded_volume_restraint()
    
    # Read stat files
    AT.read_stat_files()
    AT.write_models_info()
    AT.get_psi_stats()
    
    AT.hdbscan_clustering(["CR_sum",
                           "EV_sum",
                           "XLs_sum"])
    AT.summarize_XLs_info()
        


