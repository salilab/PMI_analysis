import IMP
import IMP.pmi
import IMP.pmi.analysis
import sys
import os

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from contact_maps import CMTable

#####################################################
# Calculate contact frequency
#####################################################
CM = CMTable(analys_dir = 'analys/',
             rmf_A = 'analys/A_models_cluster0.rmf3',
             rmf_B = 'analys/B_models_cluster0.rmf3',
             clustering_dir = 'analys/clustering/',
             cluster = 0,
             cutoff = 12.0,
             nproc = 20)

CM.add_XLs_data('XLs_data.csv',
                keys = ['Protein A',
                        'Protein B',
                        'Residue A',
                        'Residue B',
                        'Score'])
CM.compute_contact_maps()
CM.plot_contact_maps()
