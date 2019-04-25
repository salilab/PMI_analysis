import IMP
import IMP.pmi
import IMP.pmi.analysis
import sys
import os

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from contact_maps import CMTable

import glob
import numpy as np
import random

 

#####################################################
# calculate contact frequency
#####################################################
CM = CMTable(out_dir = 'analys/CMs_cluster0',
             GSMs_dir = 'analys/GSMs_2/',
             clustering_dir = 'analys/clustering_cl1/',
             cluster = 0,
             number_of_models = 15000,    
             cutoff = 10.0,
             nproc = 20)
CM.compute_contact_maps()

CM.get_close_contacts()
CM.plot_contact_maps()

