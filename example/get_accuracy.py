from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import random
import numpy as np
import pandas as pd
import sys
import multiprocessing as mp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from accuracy import *

    
selection_dictionary={"rpb1_1":[(1,846,"rpb1")],
                      "rpb1_2":[(1060,1379,"rpb1")],
		      "rpb2":[(1,1110,"rpb2")],
                      "all_sel":[(1,846,"rpb1"),(1060,1379,"rpb1"),(1,1110,"rpb2")]}


nproc = 5
refrmf = 'run_1/all_ini.rmf3'
clustering_dir = sys.argv[1]

AccuracyModels(selection_dictionary,
               clustering_dir=clustering_dir,
               ref_rmf3=refrmf,
               scores_sample_A='analys/selected_models_A_cluster0_detailed.csv',
               scores_sample_B='analys/selected_models_B_cluster0_detailed.csv',
               dir_name='run_',
               nproc=nproc)




