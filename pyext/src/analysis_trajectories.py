#!/usr/bin/env python

'''
Tools to analyze PMI MC runs
'''

from __future__ import division

import sys
import os
import math
import glob
import random
import itertools
import pandas as pd
import numpy as np
import multiprocessing as mp
from scipy import stats
from equilibration import *

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})

import seaborn as sns
import hdbscan

class AnalysisTrajectories(object):
    def __init__(self,
                 out_dirs,
                 dir_name = 'run_', 
                 analysis_dir = 'analys/',
                 nproc=6):

        self.out_dirs = out_dirs
        self.dir_name = dir_name
        self.analysis_dir = analysis_dir
        self.nproc = nproc
        self.restraint_names = {}
        self.all_fields = []

        self.th = 500
        
        # For multiprocessing
        self.manager = mp.Manager()
        self.S_all = self.manager.dict()
        self.S_files = self.manager.dict()
        self.S_dist_all = self.manager.dict()
        self.XLs_nuis = self.manager.dict()
        
        # Sample stat file
        stat_files = np.sort(glob.glob(self.out_dirs[0]+'stat.*.out'))
        self.stat2_dict = self.get_keys(stat_files[0])

        # Define with restraints to analyze
        self.Connectivity_restraint = False
        self.Excluded_volume_restraint = False
        self.XLs_restraint = False
        self.XLs_restraint_nuisances = False
        self.Multiple_XLs_restraints = False
        self.atomic_XLs_restraint = False
        self.atomic_XLs_restraint_nuisances = False
        self.Multiple_atomic_XLs_restraints = False
        
        self.EM_restraint = False
        self.Distance_restraint = False
        self.Binding_restraint = False
        self.Occams_restraint = False
        self.Occams_restraint_nuisances = False
        self.Occams_positional_restraint = False
        self.Occams_positional_nuisances = False
        self.pEMAP_restraint = False
        self.DOPE_restraint = False

        # By default, add restraints of same type
        self.sum_Connectivity_restraint = True
        self.sum_Excluded_volume_restraint = True
        self.sum_Binding_restraint = True
        self.sum_Distance_restraint =  True
        self.sum_XLs_restraint = True
        self.sum_atomic_XLs_restraint = True
        self.sum_DOPE_restraint = True
        
        self.select_XLs_satisfaction =  False
        self.select_EM_score = False
        self.select_Total_score = False
        self.Multiple_psi_values = False

        # Separate trajectories into two halves
        self.dir_halfA = np.sort(self.out_dirs)[::2]
        self.dir_halfB = np.sort(self.out_dirs)[1::2]

    def set_analyze_XLs_restraint(self,
                                  get_nuisances = True,
                                  Multiple_XLs_restraints =  False,
                                  XLs_cutoffs = {'DSSO':30.0}):
        self.XLs_restraint = True
        self.select_XLs_satisfaction = True
        if get_nuisances:
            self.XLs_restraint_nuisances = True
        if Multiple_XLs_restraints:
            self.Multiple_XLs_restraints = True
            self.sum_XLs_restraint = False
        self.XLs_cutoffs = XLs_cutoffs

    def set_analyze_atomic_XLs_restraint(self,
                                         get_nuisances = True,
                                         Multiple_atomic_XLs_restraints =  False,
                                         atomic_XLs_cutoffs = {'DSSO':30.0}):
        self.atomic_XLs_restraint = True
        self.select_atomic_XLs_satisfaction = True
        if get_nuisances:
            self.atomic_XLs_restraint_nuisances = True
        if Multiple_atomic_XLs_restraints:
            self.Multiple_atomic_XLs_restraints = True
            self.sum_atomic_XLs_restraint = False
        self.atomic_XLs_cutoffs = atomic_XLs_cutoffs
    
    def set_analyze_Connectivity_restraint(self):
        self.Connectivity_restraint = True    
       
    def set_analyze_Excluded_volume_restraint(self):
        self.Excluded_volume_restraint = True
 
    def set_analyze_EM_restraint(self):
        self.EM_restraint = True

    def set_analyze_Distance_restraint(self):
        self.Distance_restraint = True

    def set_analyze_Binding_restraint(self):
        self.Binding_restraint = True
        
    def set_analyze_Occams_restraint(self):
         self.Occams_restraint = True
         
    def set_analyze_Occams_positional_restraint(self):
        self.Occams_positional_restraint = True
        self.Occams_positional_nuisances = True
    
    def set_analyze_pEMAP_restraint(self):
        self.pEMAP_restraint = True

    def set_analyze_DOPE_restraint(self):
        self.DOPE_restraint = True

    def set_select_by_XLs_satisfation(self):
        self.select_EM_score = True
        
    def set_select_by_Total_score(self, score_cutoff):
        self.select_Total_score = True
        self.cutoff_Total_score = score_cutoff

    def set_select_by_EM_score(self, score_cutoff):
        self.select_EM_score = True
        self.cutoff_EM_score = score_cutoff
    
    def get_restraint_fields(self):
        '''
        Gather information to retrieve from
        stat file for each trajectory
        '''
    
        # Get the field number for fframe and total score
        self.rmf_file_field = self.get_field_id(self.stat2_dict, 'rmf_file')
        rmf_frame_index = self.get_field_id(self.stat2_dict, 'rmf_frame_index')
        fframe = self.get_field_id(self.stat2_dict, 'MonteCarlo_Nframe')
        Total_score = self.get_field_id(self.stat2_dict, 'Total_Score')
    
        # Restraint names
        self.all_fields = fframe+rmf_frame_index+Total_score
        self.restraint_names[0] = 'MC_frame'
        self.restraint_names[1] = 'rmf_frame_index'
        self.restraint_names[2] = 'Total_score'
        name_i = 3

        # Percent satisfaction fields
        self.percent_satisfied = []
        
        # All possible restraints
        if self.Excluded_volume_restraint:
            EV = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('ExcludedVolumeSphere_' in self.stat2_dict[k])}
            if len(EV) == 1:
                self.restraint_names[name_i] = 'EV_sum'
                name_i += 1
            else:
                for k, v in EV.items():
                    self.restraint_names[name_i] = 'EV_'+k.split('ExcludedVolumeSphere_')[-1]
                    name_i += 1
            self.all_fields += EV.values()

        if self.Connectivity_restraint:
            CR = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('ConnectivityRestraint_' in self.stat2_dict[k])}
            if len(CR) == 1:
                self.restraint_names[name_i] = 'CR_sum'
                name_i += 1
            else:
                for k, v in CR.items():
                    self.restraint_names[name_i] = 'CR_'+k.split('ConnectivityRestraint_')[-1]
                    name_i += 1
            self.all_fields += CR.values()  

        if self.XLs_restraint:
            XLs = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('CrossLinkingMassSpectrometryRestraint_Data_Score' in self.stat2_dict[k])}
            if len(XLs) == 1:
                self.restraint_names[name_i] = 'XLs_sum'
                name_i += 1
            else:
                for k, v in XLs.items():
                    self.restraint_names[name_i] = 'XLs_'+k.split('CrossLinkingMassSpectrometryRestraint_')[-1]
                    name_i += 1
            self.all_fields += XLs.values()

            # Get Psi score
            Psis = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('CrossLinkingMassSpectrometryRestraint_PriorPsi_Score' in self.stat2_dict[k])}
            if len(Psis) == 1:
                self.restraint_names[name_i] = 'Psi_scores_Sum'
                name_i += 1
            else:
                for k, v in Psis.items():
                    self.restraint_names[name_i] = 'Psi_scores_'+k.split('CrossLinkingMassSpectrometryRestraint_PriorPsi_Score_')[-1]
                    name_i += 1
            self.all_fields += Psis.values()
            
            # Other quantities associated to XLs (distances and nuisances)
            XLs_dist = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('CrossLinkingMassSpectrometryRestraint_Distance_' in self.stat2_dict[k])}
            self.XLs_info = XLs_dist
            if self.XLs_restraint_nuisances:
                XLs_nuis = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('CrossLinkingMassSpectrometryRestraint_Psi_' in self.stat2_dict[k] and 'MonteCarlo_' not in self.stat2_dict[k])}
                self.psi_head = self.get_str_match(list(XLs_nuis.keys()))
                for k, v in XLs_nuis.items():
                    self.restraint_names[name_i] = 'Psi_vals_'+k.split('CrossLinkingMassSpectrometryRestraint_Psi_')[-1]
                    name_i += 1
                self.all_fields += XLs_nuis.values()

                if len(XLs_nuis.keys())> 1:
                    self.Multiple_psi_values = True
                self.XLs_info.update(XLs_nuis)
                mm = sorted(list([v for v in XLs_nuis.keys()]))
                ss = sorted(list([v+'_std' for v in XLs_nuis.keys()]))
                self.DF_XLs_psi = pd.DataFrame(columns=['Trajectory']+mm+ss)
            self.xls_fields = self.XLs_info.values()

        # Atomic XLs restraint
        if self.atomic_XLs_restraint:
            atomic_XLs = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('AtomicXLRestraint_Score' in self.stat2_dict[k])}
            if len(atomic_XLs) == 1:
                self.restraint_names[name_i] = 'atomic_XLs_sum'
                name_i += 1
            else:
                for k, v in XLs.items():
                    self.restraint_names[name_i] = 'atomic_XLs_'+k.split('AtomicXLRestraint_Score_')[-1]
                    name_i += 1
            self.all_fields += atomic_XLs.values()

            # Get Psi score
            Psi_score = self.get_field_id(self.stat2_dict,'AtomicXLRestraint_psi_Score')
            self.restraint_names[name_i] = 'atomic_Psi_score'
            name_i += 1
            self.all_fields += Psi_score
            
            # Other quantities associated to XLs (distances and nuisances)
            atomic_XLs_dist = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('AtomicXLRestraint_' in self.stat2_dict[k] and 'BestDist' in self.stat2_dict[k])}
            self.atomic_XLs_info = atomic_XLs_dist
            if self.atomic_XLs_restraint_nuisances:
                atomic_XLs_nuis = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('AtomicXLRestraint_psi_Score' in self.stat2_dict[k] and 'MonteCarlo_' not in self.stat2_dict[k])}
                if len(atomic_XLs_nuis.keys())> 1:
                   self.Multiple_atomic_psi_values = True
                self.atomic_XLs_info.update(atomic_XLs_nuis)
                mm = sorted(list([v for v in atomic_XLs_nuis.keys()]))
                ss = sorted(list([v+'_std' for v in atomic_XLs_nuis.keys()]))
                self.DF_atomic_XLs_psi = pd.DataFrame(columns=['Trajectory']+mm+ss)
            self.xls_fields = self.atomic_XLs_info.values()
            
        if self.EM_restraint:
            EM3D = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('GaussianEMRestraint_' in self.stat2_dict[k] and 'sigma' not in self.stat2_dict[k])}
            if len(EM3D) == 1:
                self.restraint_names[name_i] = 'EM3D_sum'
                name_i += 1
            else:
                for k, v in EM3D.items():
                    self.restraint_names[name_i] = 'EM3D_'+k.split('GaussianEMRestraint_')[-1]
                    name_i += 1
            self.all_fields += EM3D.values()

        if self.Occams_restraint:
            Occ = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('OccamsRestraint_0' in self.stat2_dict[k] and 'sigma' not in self.stat2_dict[k])}
            if len(Occ) == 1:
                self.restraint_names[name_i] = 'Occams_sum'
                name_i += 1
            else:
                for k, v in Occ.items():
                    self.restraint_names[name_i] = 'Occams_'+k.split('OccamsRestraint_')[-1]
                    name_i += 1
            self.all_fields += Occ.values()

            # Percent of restraints satisfied
            self.Occams_satif = self.get_field_id(self.stat2_dict, 'OccamsRestraint_satisfied_0')
            
        if self.pEMAP_restraint:
            pEMAP = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('SimplifiedPEMAP_Score_' in self.stat2_dict[k] and ':' not in self.stat2_dict[k])}
            if len(pEMAP) == 1:
                self.restraint_names[name_i] = 'pEMAP_sum'
                name_i += 1
            else:
               for k, v in pEMAP.items():
                self.restraint_names[name_i] = 'pEMap_'+k.split('_')[-1]
                name_i += 1
            self.all_fields += pEMAP.values()
            
            # Percent of restraints satisfied
            self.pEMAP_satif = self.get_field_id(self.stat2_dict, 'SimplifiedPEMAP_Satisfied')
            # All pE-MAP distances
            pEMAP_dist = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('SimplifiedPEMAP_Distance_' in self.stat2_dict[k])}
            
            # Get target distance from stat file
            d_pEMAP = []
            v_pEMAP = []
            for k, val in pEMAP_dist.items():
                d_pEMAP.append(float(k.split('_')[-1]))        
                v_pEMAP.append(val)

        if self.Distance_restraint:
            DR = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('DistanceRestraint_Score_' in self.stat2_dict[k])}
            if len(DR) == 1:
                self.restraint_names[name_i] = 'DR_sum'
                name_i += 1
            else:
                for k, v in DR.items():
                    self.restraint_names[name_i] = 'DR_'+k.split('DistanceRestraint_Score_')[-1]
                    name_i += 1
            self.all_fields += DR.values()

        if self.Binding_restraint:
            BR = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('ResidueBindingRestraint_score_' in self.stat2_dict[k])}
            if len(BR) == 1:
                self.restraint_names[name_i] = 'BR_sum'
                name_i += 1
            else:
                for k, v in BR.items():
                    self.restraint_names[name_i] = 'BR_'+k.split('ResidueBindingRestraint_score_')[-1]
                    name_i += 1
            self.all_fields += BR.values()

        if self.DOPE_restraint:
            DOPE = {self.stat2_dict[k]: k for k in self.stat2_dict.keys() if ('DOPE_Restraint_score' in self.stat2_dict[k])}
            if len(DOPE) == 1:
                self.restraint_names[name_i] = 'DOPE_sum'
                name_i += 1
            else:
                for k, v in DOPE.items():
                    self.restraint_names[name_i] = 'DOPE_'+k.split('DOPERestraint_score_')[-1]
                    name_i += 1
            self.all_fields += DOPE.values()

    def read_DB(self, db_file):
        ''' Read database '''
        DB = {}
        i = 0
        for line in open(db_file):
            vals =  line.split()
            DB[str(i)] = vals
            i  += 1
        return DB

    def get_keys(self, stat_file):
        ''' Get all keys in stat file '''
        for line in open(stat_file).readlines():
            d = eval(line)
            klist = list(d.keys())
            # check if it is a stat2 file
            if "STAT2HEADER" in klist:
                import operator
                isstat2 = True
                for k in klist:
                    if "STAT2HEADER" in str(k):
                        del d[k]
                stat2_dict = d
                # get the list of keys sorted by value
                kkeys = [k[0]
                         for k in sorted(stat2_dict.items(), key=operator.itemgetter(1))]
                klist = [k[1]
                         for k in sorted(stat2_dict.items(), key=operator.itemgetter(1))]
                invstat2_dict = {}
                for k in kkeys:
                    invstat2_dict.update({stat2_dict[k]: k})
            else:
                isstat1 = True
                klist.sort()
            break
        
        return stat2_dict

    def read_stat_files(self):
        
        # Split directories
        ND = int(np.ceil(len(self.out_dirs)/float(self.nproc)))
        out_dirs_dict = {}
        for k in range(self.nproc-1):
            out_dirs_dict[k] = list(self.out_dirs[(k*ND):(k*ND+ND)])
        out_dirs_dict[self.nproc-1] = list(self.out_dirs[((self.nproc-1)*ND):(len(self.out_dirs))])

        # Define an output queue
        output = mp.Queue()
        
        # Setup a list of processes that we want to run
        processes = [mp.Process(target=self.read_traj_info, args=((out_dirs_dict[x],))) for x in range(self.nproc)]

        # Run processes
        for p in processes:
            p.start()
            
        # Exit the completed processes
        for p in processes:
            p.join()

    def read_stats_detailed(self, traj, stat_files, query_fields, dist_fields, satif_fields, query_rmf_file):
        """
        Detailed reading of stats files that includes
        the rmf in which the frame is.
        To be used when using rmf_slice
        """
        
        S_scores = []
        S_dist = []
        P_satif = []
        #frames_dic = {}
        for file in stat_files:
            line_number=0
            for line in open(file).readlines():
                line_number += 1
                try:
                    d = eval(line)
                except:
                    print("# Warning: skipped line number " + str(line_number) + " not a valid line")
                    break 

                if line_number > 1:
                    frmf = [d[field] for field in query_rmf_file][0]
                    s0 = [float(d[field]) for field in self.all_fields]+[traj, frmf]
                    S_scores.append(s0)
                    if self.XLs_restraint or  self.atomic_XLs_restraint:
                        d0 = [s0[0]] + [float(d[field]) for field in self.xls_fields]
                        S_dist.append(d0)
                    if satif_fields:
                        p0 = [s0[0]] + [float(d[field]) for field in satif_fields]
                        P_satif.append(p0)
                    
        # Sort based on frame
        S_scores.sort(key=lambda x: float(x[0]))
        
        # Convert into pandas DF
        column_names = [self.restraint_names[x] for x in sorted(self.restraint_names.keys())] + ['traj', 'rmf3_file']
        DF = pd.DataFrame(S_scores, columns = column_names)
        
        # If some restraints need to be added
        if self.XLs_restraint and self.sum_XLs_restraint:
            XLs_sum = pd.Series(self.add_restraint_type(DF, 'XLs_'))
            DF = DF.assign(XLs_sum=XLs_sum.values)

        if self.atomic_XLs_restraint and self.sum_atomic_XLs_restraint:
            atomic_XLs_sum = pd.Series(self.add_restraint_type(DF, 'atomic_XLs_'))
            DF = DF.assign(atomic_XLs_sum=atomic_XLs_sum.values)
       
        if self.Excluded_volume_restraint and self.sum_Excluded_volume_restraint:
            EV_sum = pd.Series(self.add_restraint_type(DF, 'EV_'))
            DF = DF.assign(EV_sum=EV_sum.values)
 
        if self.Connectivity_restraint and self.sum_Connectivity_restraint:
            CR_sum = pd.Series(self.add_restraint_type(DF, 'CR_'))
            DF = DF.assign(CR_sum=CR_sum.values)
            
        if self.Binding_restraint and self.sum_Binding_restraint:
            BR_sum = pd.Series(self.add_restraint_type(DF, 'BR_'))
            DF = DF.assign(BR_sum=BR_sum.values)
            
        if self.Distance_restraint and self.sum_Distance_restraint:
            DR_sum = pd.Series(self.add_restraint_type(DF, 'DR_'))
            DF = DF.assign(DR_sum=DR_sum.values)
            
        # Get distance fields
        if self.XLs_restraint:
            S_dist = np.array(S_dist)
            S_dist = S_dist[S_dist[:,0].argsort()]
            # Convert in DF
            column_names_XLs = ['MC_frame'] + list(self.XLs_info.keys())
            DF_XLs = pd.DataFrame(S_dist, columns = column_names_XLs)
        if self.atomic_XLs_restraint:
            S_dist = np.array(S_dist)
            S_dist = S_dist[S_dist[:,0].argsort()]
            # Convert in DF
            column_names_XLs = ['MC_frame'] + list(self.atomic_XLs_info.keys())
            DF_atomic_XLs = pd.DataFrame(S_dist, columns = column_names_XLs)
     
        if satif_fields:  
            P_satif = np.array(P_satif)
            P_satif = P_satif[P_satif[:,0].argsort()]

        #return DF, DF_XLs, P_satif, frames_dic
        if self.XLs_restraint:
            return DF, DF_XLs, P_satif
        elif self.atomic_XLs_restraint:
            return DF, DF_atomic_XLs, P_satif
        elif satif_fields:
            return DF, None, P_satif
        else:
            return DF, None, None

    def add_restraint_type(self, DF, key_id):
        temp_fields = [v for v in DF.columns.values if key_id in v]
        DF_t = DF[temp_fields]
        DF_s = DF_t.sum(axis=1)
        return DF_s
    
    def read_traj_info(self, out_dirs_sel):
        # Dictionary to put the scores of all the trajectories
        # files_dic: keys: frames, values: rmf3 file

        #rmf_file = get_field_id(stat2_dict, 'rmf_file')
        #rmf_frame_index = get_field_id(stat2_dict, 'rmf_frame_index')
        if isinstance(out_dirs_sel, str):
            out_dirs_sel = [out_dirs_sel]
        
        for out in out_dirs_sel:
            traj = [x for x in out.split('/') if self.dir_name in x][0]
            traj_number = int(traj.split(self.dir_name)[1])
            stat_files = np.sort(glob.glob(out+'stat.*.out'))
            #if self.XLs_restraint:
            #    S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
            #                                                             stat_files,
            #                                                             self.all_fields,
            #                                                             self.XLs_info.values(),
            #                                                             None,
            #                                                             self.rmf_file_field)
                
                
            if self.atomic_XLs_restraint and not self.pEMAP_restraint and not self.Occams_restraint:
                S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
                                                                         stat_files,
                                                                         self.all_fields,
                                                                         self.atomic_XLs_info.values(),
                                                                         None,
                                                                         self.rmf_file_field)
                
            elif self.pEMAP_restraint and not self.XLs_restraint:
                S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
                                                                         stat_files,
                                                                         self.all_fields,
                                                                         None,
                                                                         self.pEMAP_satif,
                                                                         self.rmf_file_field)
    
            elif self.pEMAP_restraint and self.XLs_restraint:
                S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
                                                                         stat_files,
                                                                         self.all_fields,
                                                                         self.XLs_info.values(),
                                                                         self.pEMAP_satif,
                                                                         self.rmf_file_field)
            elif self.Occams_restraint  and not self.XLs_restraint:
                S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
                                                                         stat_files,
                                                                         self.all_fields,
                                                                         None,
                                                                         self.Occams_satif,
                                                                         self.rmf_file_field)
    
            elif self.Occams_restraint and self.XLs_restraint:
                S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
                                                                         stat_files,
                                                                         self.all_fields,
                                                                         self.XLs_info.values(),
                                                                         self.Occams_satif,
                                                                         self.rmf_file_field)

            
                
            else:
                S_tot_scores, S_dist, P_satif = self.read_stats_detailed(traj,
                                                                         stat_files,
                                                                         self.all_fields,
                                                                         None,
                                                                         None,
                                                                         self.rmf_file_field)


            print('The mean score, min score, and n frames are: ', traj_number,
                  np.mean(S_tot_scores['Total_score'].iloc[self.th:]),
                  np.min(S_tot_scores['Total_score'].iloc[self.th:]),
                  len(S_tot_scores))
                
            # Selection of just sums (default)
            sel_entries = ['Total_score']+[v for v in S_tot_scores.columns.values if 'sum' in v]
            
            # If specified by user, can look at individual contributions
            if self.Connectivity_restraint== True and self.sum_Connectivity_restraint == False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'CR_' in v and 'sum' not in v]
            if self.Excluded_volume_restraint== True and self.sum_Excluded_volume_restraint == False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'EV_' in v and 'sum' not in v]
            if self.Binding_restraint == True and self.sum_Binding_restraint == False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'BR_' in v and 'sum' not in v]
            if self.Distance_restraint ==  True and self.sum_Distance_restraint ==  False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'DR_' in v and 'sum' not in v]
            if self.XLs_restraint == True and self.sum_XLs_restraint == False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'XLs_' in v and 'sum' not in v]
            if self.atomic_XLs_restraint == True and self.sum_atomic_XLs_restraint == False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'atomic_XLs_' in v and 'sum' not in v]
            if self.DOPE_restraint == True and self.sum_DOPE_restraint == False:
                sel_entries += [v for v in S_tot_scores.columns.values if 'DOPE_' in v and 'sum' not in v]

            # Also add nuisances parameter
            sel_entries += [v for v in S_tot_scores.columns.values if 'Psi' in v and 'sum' not in v]
            
            # Detect equilibration time
            ts_eq = []
            for r in sel_entries:
                try:
                    [t, g, N] = detectEquilibration(np.array(S_tot_scores[r].loc[self.th:]), nskip=5, method='multiscale')
                    ts_eq.append(t)
                except:
                    ts_eq.append(0)
            print('ts_eq', ts_eq)
            ts_max = np.max(ts_eq)+self.th
        
            # Plot the scores and restraint satisfaction
            file_out = 'plot_scores_%s.pdf'%(traj_number)               
            self.plot_scores_restraints(S_tot_scores[['MC_frame']+sel_entries], ts_eq, file_out)
        
            if self.pEMAP_restraint:
                file_out_pemap = 'plot_pEMAP_%s.pdf'%(traj_number) 
                self.plot_pEMAP_distances(P_satif, file_out_pemap)

            if self.Occams_restraint:
                file_out_occams = 'plot_Occams_satisfaction_%s.pdf'%(traj_number) 
                self.plot_Occams_satisfaction(P_satif, file_out_occams)
                
            # Check how many XLs are satisfied
            if self.XLs_restraint:
                S_tot_scores, S_dist = self.analyze_XLs_values(S_tot_scores,
                                                               S_dist,
                                                               self.XLs_cutoffs,
                                                               atomic_XLs = False,
                                                               traj_number = traj_number,
                                                               ts_max = ts_max )
            
            if self.atomic_XLs_restraint:
                S_tot_scores, S_dist = self.analyze_XLs_values(S_tot_scores,
                                                               S_dist,
                                                               self.atomic_XLs_cutoffs,
                                                               atomic_XLs = True,
                                                               traj_number = traj_number,
                                                               ts_max = ts_max )
                                                                                                 
            # Add scores to dictionary XLs_satif
            
            # Add half info
            if out in self.dir_halfA:
                S_tot_scores = S_tot_scores.assign(half = pd.Series(['A']*len(S_tot_scores), index=S_tot_scores.index).values)
            elif out in self.dir_halfB:
                S_tot_scores = S_tot_scores.assign(half = pd.Series(['B']*len(S_tot_scores), index=S_tot_scores.index).values)
            else:
                S_tot_scores = S_tot_scores.assign(half = pd.Series([0]*len(S_tot_scores), index=S_tot_scores.index).values)

            self.S_all[out] = S_tot_scores[ts_max:]
            if self.XLs_restraint:
                self.S_dist_all[out] = S_dist[ts_max:]
            if self.atomic_XLs_restraint:
                self.S_dist_all[out] = S_dist[ts_max:]
                
    def get_field_id(self, dict,val):
        ''' 
        For single field, get number of fields in stat file 
        '''
        return [k for k in dict.keys() if dict[k]==val]

    def get_XLs_satisfaction(self, S_dist, XLs_cutoffs, atomic_XLs, type_XLs = None, type_psi = None):
        if type_XLs and not type_psi:
            dist_columns = [x for x in S_dist.columns.values if ('Distance_' in x and type_XLs in x)]
            cutoff = XLs_cutoffs[type_XLs]
        elif type_psi and not type_XLs:
            dist_columns = [x for x in S_dist.columns.values if ('Distance_' in x and type_psi in x)]
            cutoff = list(XLs_cutoffs.values())[0]
        elif type_XLs and type_psi:
            dist_columns = [x for x in S_dist.columns.values if ('Distance_' in x and type_XLs in x and type_psi in x)]
            cutoff = XLs_cutoffs[type_XLs]
        elif atomic_XLs:
            dist_columns = [x for x in S_dist.columns.values if ('BestDist' in x)]
            cutoff = list(XLs_cutoffs.values())[0]
        else:
            dist_columns = [x for x in S_dist.columns.values if 'Distance_' in x]
            cutoff = list(XLs_cutoffs.values())[0]
        XLs_dist = np.array(S_dist[dist_columns])
        
        perc_per_step = []
        for row in XLs_dist:
            n_satif = 0
            for i in row:
                if i <= cutoff:
                    n_satif += 1
            perc_per_step.append(float(n_satif)/len(row))
        return perc_per_step
    
    def plot_scores_restraints(self, selected_scores, ts_eq, file_out):
        '''
        For each trajectory plot all restraint scores
        '''
        n_bins=20
        ts_max = np.max(ts_eq)
        n_res = len(selected_scores.columns.values)-1
    
        fig, ax = pl.subplots(figsize=(2.0*n_res, 4.0), nrows=2, ncols=n_res)
        axes = ax.flatten()
        for i, c in enumerate(selected_scores.columns.values[1:]):
            axes[i].plot(selected_scores['MC_frame'].loc[self.th::10], selected_scores[c].loc[self.th::10], color='b',alpha=0.5)
            axes[i].axvline(ts_eq[i], color='grey')
            axes[i].set_title(c, fontsize=14)
            axes[i].set_xlabel('Step',fontsize=12)
            if i == 0:
                axes[i].set_ylabel('Score (a.u.)',fontsize=12)
        
        for i, c in enumerate(selected_scores.columns.values[1:]):
            axes[i+n_res].hist(selected_scores[c].loc[ts_eq[i]::10], n_bins, histtype='step',fill=False, color='orangered',alpha=0.9)
            axes[i+n_res].hist(selected_scores[c].loc[ts_max::10], n_bins, histtype='step',fill=False, color='gold',alpha=0.9)
            axes[i+n_res].set_xlabel('Score (a.u.)',fontsize=12)
            if i == 0:
                axes[i+n_res].set_ylabel('Density',fontsize=12)

        pl.tight_layout(pad=0.5, w_pad=0.1, h_pad=2.0)
        fig.savefig(self.analysis_dir+file_out) 

    def analyze_XLs_values(self, S_tot_scores, S_dist, XLs_cutoffs, atomic_XLs, traj_number, ts_max):

        sel_nuis = [v for v in S_dist.columns.values if 'Psi' in v]
                
        # Get XLS satisfaction, append to S_dist
        XLs_satif_fields = []
        if self.Multiple_XLs_restraints:
            for type_XLs in self.XLs_cutoffs.keys():
                XLs_satif = self.get_XLs_satisfaction(S_dist, XLs_cutoffs, atomic_XLs, type_XLs = type_XLs)
                temp_name = 'XLs_satif_'+type_XLs.rstrip()
                S_tot_scores = S_tot_scores.assign(XLs_satif=pd.Series(XLs_satif))
                S_tot_scores.rename(columns={'XLs_satif':temp_name}, inplace=True)
                XLs_satif_fields.append(temp_name)
        elif self.Multiple_psi_values:
            all_psis = [v.split(self.psi_head)[1] for v in  self.DF_XLs_psi.columns.values[1:] if 'std' not in v]
            for type_psi in all_psis:
                XLs_satif = self.get_XLs_satisfaction(S_dist, XLs_cutoffs, atomic_XLs, type_psi = type_psi)
                temp_name = 'XLs_satif_'+type_psi
                S_tot_scores = S_tot_scores.assign(XLs_satif=pd.Series(XLs_satif))
                S_tot_scores.rename(columns={'XLs_satif':temp_name}, inplace=True)
                XLs_satif_fields.append(temp_name)
        else:
            XLs_satif = self.get_XLs_satisfaction(S_dist, XLs_cutoffs, atomic_XLs)
            S_tot_scores = S_tot_scores.assign(XLs_satif=pd.Series(XLs_satif))
            XLs_satif_fields.append('XLs_satif')
            
        file_out_xls = 'plot_XLs_%s.pdf'%(traj_number)
        self.plot_XLs_satisfaction(S_tot_scores['MC_frame'].values, S_tot_scores[XLs_satif_fields], S_dist[sel_nuis], ts_max, file_out_xls)
        
        # Add percent of satisfied XLs to DF
        if self.XLs_restraint_nuisances:
            nuis_mean = []
            nuis_std = []
            nuis_fields = sorted([n for n in S_dist.columns.values if 'Psi' in n])
            for nuis in sorted(nuis_fields):
                # Create df instead of dic
                nuis_mean.append(np.mean(S_dist[nuis].loc[ts_max:]))
                nuis_std.append(np.std(S_dist[nuis].loc[ts_max:]))
            self.XLs_nuis[traj_number] = list(nuis_mean) + list(nuis_std)

        return S_tot_scores, S_dist
    
    def plot_XLs_satisfaction(self, t, perc_per_step, nuis_vals, ts_max, file_out):
        c = ['gold', 'orange', 'red', 'blue', 'green']
        n_bins = 20
        
        fig, ax = pl.subplots(figsize=(10.0, 4.0), nrows=1, ncols=3)
        axes = ax.flatten()
        for i, c in enumerate(perc_per_step.columns.values):
            label = c
            axes[0].plot(t[::100], perc_per_step[c].loc[::100], label=label)
        axes[0].set_title('XLs restraint satisfaction', fontsize=14)
        axes[0].set_xlabel('Step',fontsize=12)
        axes[0].set_ylabel('Percent Satisfied',fontsize=12)
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(handles[::-1], labels[::-1])

        for i, c in enumerate(nuis_vals.columns.values):
            label = c.split('CrossLinkingMassSpectrometryRestraint_')[-1]
            axes[1].plot(t[::100], nuis_vals[c].loc[::100], label=label)
        axes[1].set_title('Psi nuisance parameters', fontsize=14)
        axes[1].set_xlabel('Step',fontsize=12)
        axes[1].set_ylabel('Psi',fontsize=12)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[1].legend(handles[::-1], labels[::-1])

        for i,c in enumerate(nuis_vals.columns.values):
            label = c.split('CrossLinkingMassSpectrometryRestraint_')[-1]
            axes[2].hist(nuis_vals[c].loc[ts_max:],n_bins, histtype='step',fill=False, label=label)
            
        axes[2].set_title('Psi nuisance parameters', fontsize=14)
        axes[2].set_xlabel('Psi',fontsize=12)
        axes[2].set_ylabel('Density',fontsize=12)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[2].legend(handles[::-1], labels[::-1])
        
        pl.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.5)
        fig.savefig(self.analysis_dir+file_out)

    def plot_selected_fields(self, field1, field2, k=None):
        from sklearn.cluster import KMeans
        from sklearn.externals import joblib

        s1 = [self.S_all[dd][field1] for dd in self.S_all.keys()]
        s2 = [self.S_all[dd][field2] for dd in self.S_all.keys()]

        # Now, for GSMs
        traj_0 = self.S_all.keys()[0]
        col_sel = [v for v in self.S_all[traj_0].columns.values if field1 in v or field2 in v]
        i=0
        for traj in self.S_all.keys():
            T = self.S_all[traj]
            # Get MC frames for GSMs A and B
            frames_sel = self.gsms_frames[traj]
            if i == 0 and len(frames_sel) > 0:
                all_GSMs_sel = T.loc[T['MC_frame'].isin(frames_sel)][col_sel]
                i =  1
            elif i==1 and len(frames_sel) > 0:
                all_GSMs_sel.append(T.loc[T['MC_frame'].isin(frames_sel)][col_sel], ignore_index = True)
        all_GSMs_sel = np.array(all_GSMs_sel)
        
        d1 =  pd.concat(s1)
        d2 =  pd.concat(s2)

        D = np.array([d1,d2]).T
        
        file_out = 'scores_'+field1+'_'+field2+'.pdf'
        xi = np.min(D[:,0])
        xe = np.max(D[:,0])
        yi = np.min(D[:,1])
        ye = np.max(D[:,1])

        fig, ax = pl.subplots(figsize=(10.0, 4.0), nrows=1, ncols=3)
        ax[0].scatter(D[::20,0], D[::20,1],alpha=0.8, s=3)
        ax[0].set_title('Scores distributions', fontsize=14)
        ax[0].set_xlabel(field1+' (a.u.)',fontsize=12)
        ax[0].set_ylabel(field2+' (a.u.)',fontsize=12)
        ax[0].set_xlim(xi,xe)
        ax[0].set_ylim(yi,ye)

        # Kmeans
        cmap='viridis'
        kmeans = KMeans(n_clusters=5)
        kmeans.fit(D)
        y_kmeans = kmeans.predict(D[::20])
        centers = kmeans.cluster_centers_
       
        ax[1].scatter(D[::20,0],D[::20,1], c=y_kmeans, s=3, cmap=cmap)
        ax[1].scatter(centers[:,0], centers[:,1],c='grey',s=100, alpha=0.75)
        ax[1].set_title('Clustered scores distributions (KMeans)', fontsize=14)
        ax[1].set_xlabel(field1+' (a.u.)',fontsize=12)
        ax[1].set_ylabel(field2+' (a.u.)',fontsize=12)
        ax[1].set_xlim(xi,xe)
        ax[1].set_ylim(yi,ye)
        
        ax[2].scatter(all_GSMs_sel[::20,0],all_GSMs_sel[::20,1],alpha=0.8, s=3)
        ax[2].set_title('Scores distributions', fontsize=14)
        ax[2].set_xlabel(field1+' (a.u.)',fontsize=12)
        ax[2].set_ylabel(field2+' (a.u.)',fontsize=12)
        ax[2].set_xlim(xi,xe)
        ax[2].set_ylim(yi,ye)
        
        pl.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.5)
        fig.savefig(self.analysis_dir+file_out)

        # Save model
        filename = self.analysis_dir+'/finalized_model_'+field1+'_'+field2+'.dat'
        joblib.dump(kmeans, filename)

        print('Shape D, GSMs', np.shape(D), np.shape(all_GSMs_sel))
        
    def get_Psi_stats(self, atomic_XLs = False):
        '''
        Organize Psi values into DataFrame
        Get mean value for extracting models 
        (Does not work for atomic XLs restraint)
        '''
        if atomic_XLs:
            DF_XLs_psi  =  self.DF_atomic_XLs_psi
        else:
            DF_XLs_psi  =  self.DF_XLs_psi

        
        for k,v in self.XLs_nuis.items():
            DF_XLs_psi = DF_XLs_psi.append(pd.Series([k]+v, index = DF_XLs_psi.columns.values), ignore_index=True) 

        DF_XLs_psi.to_csv(self.analysis_dir+'XLs_all_Psi.csv')

        psi_cols =  DF_XLs_psi.columns.values[1:]
        psi_vals = DF_XLs_psi[psi_cols]
        if self.Multiple_XLs_restraints or self.Multiple_psi_values:
            self.psi_mean = psi_vals.mean()
        else:
            self.psi_mean = psi_vals.mean().mean()
           
        print('The average XLs Psi parameter is: ', self.psi_mean)

    def write_models_info(self):
        ''' Write info of all models after equilibration'''
        
        for k, T in self.S_all.items():
            kk = k.split(self.dir_name)[-1].split('/')[0]
            T.to_csv(self.analysis_dir+'all_info_'+str(kk)+'.csv')

    def read_models_info(self):
        ''' Read info of all models after equilibration'''
        
        info_files = glob.glob(self.analysis_dir+'all_info_*.csv')
        for f in info_files:
            k = f.split('all_info_')[-1].split('.csv')[0]
            df = pd.read_csv(f)
            self.S_all[k] = df
    
    def do_GSMs_selection(self, dir_half):
        '''
        Select GSMs based on a series (or combination) of conditions.
        Not all possible conditions/combinations are implemented
        Need to extend
        '''
    
        good_half = []
        for dd in dir_half:
            T = self.S_all[dd]
            all_masks = []
            if self.select_XLs_satisfaction == True:
                
                if self.Multiple_psi_values == True:
                    print('Selecting models based on multiple Psi values')
                    XLs_psi_cutoffs = {}
                    all_psis = [v.split('_Psi_')[1] for v in  self.DF_XLs_psi.columns.values[1:] if 'std' not in v]
                    for psi_id in all_psis:
                        satif_field = [field for field in T.columns.values if (psi_id in field and 'satif' in field)][0]
                        psi_field = [field for field in self.psi_mean.index.values if psi_id in field if 'std' not in field]
                        psi_val = float(self.psi_mean[psi_field])
                        psi_std = float(self.psi_mean[psi_field[0]+'_std'])
                        XLs_psi_cutoffs[satif_field] = 1.0-7.0*psi_val
                    # Selected XLs by different psi
                    vv = [XLs_psi_cutoffs[v] for v in XLs_psi_cutoffs.keys()]
                    mask_XLs = (T[list(XLs_psi_cutoffs.keys())] >= pd.Series(vv, index=XLs_psi_cutoffs.keys())).all(axis=1)
                    all_masks.append(mask_XLs)

                # Check it still works with multiple differernt XL restraints
                #elif self.Multiple_psi_values == True and self.Multiple_XLs_restraints == False:
                #elif self.Multiple_psi_values == False and self.Multiple_XLs_restraints == True:

                else:
                   XLs_cutoffs = {}
                   XLs_cutoffs['XLs_sum'] = (1-self.psi_mean)
                   mask_XLs = (T[XLs_cutoffs.keys()] <= pd.Series(XLs_cutoffs)).all(axis=1)
                   all_masks.append(mask_XLs) 

            if self.select_EM_score == True:
                EM_cutoffs = {}
                EM_cutoffs['EM3D_sum'] = self.cutoff_EM_score
                mask_EM = (T[EM_cutoffs.keys()] <= pd.Series(EM_cutoffs)).all(axis=1)
                all_masks.append(mask_EM)
                
            if self.select_Total_score == True:
                Total_score_cutoffs = {}
                Total_score_cutoffs['Total_score'] = self.cutoff_Total_score
                mask_Total = (T[Total_score_cutoffs.keys()] <= pd.Series(Total_score_cutoffs)).all(axis=1)
                all_masks.append(mask_Total)

            # Select based on the intersection of all masks
            if len(all_masks)>0:
                mask = pd.concat(all_masks, axis = 1).all(axis=1)
                sele = T.loc[mask,:]
                good_half.append(sele)
            else:
                sele = T
                good_half.append(sele)

            # Add frames to dictionary
            if dd in self.gsms_frames.keys():
                self.gsms_frames[dd].append(sele['MC_frame'])   
            else:
                self.gsms_frames[dd] = sele['MC_frame']    

        return good_half

    def do_hdbscan_clustering(self,
                              selected_scores,
                              min_cluster_size=150,
                              min_samples=5,
                              skip=1):

        '''
        DO HDBSCAN clustering for selected restraint and/or nuisance parameters
        '''
 
        all_dfs = [self.S_all[dd] for dd in self.S_all.keys()]
        S_comb = pd.concat(all_dfs)
        
        S_comb_sel = S_comb[selected_scores].iloc[::skip]
        S_comb_all = S_comb.iloc[::skip]       
 
        hdbsc = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size,
                                min_samples=min_samples).fit(S_comb_sel)
        labels = hdbsc.labels_

        # Add clusters labels
        S_comb_sel = S_comb_sel.assign(cluster = pd.Series(hdbsc.labels_, index=S_comb_sel.index).values)
        S_comb_all = S_comb_all.assign(cluster = pd.Series(hdbsc.labels_, index=S_comb_all.index).values)
        
        print('Number of unique clusters: ', len(np.unique(hdbsc.labels_)))
        
        # Write and plot info from clustering
        self.plot_hdbscan_clustering(S_comb_sel, selected_scores)
        self.do_write_hdbscan_clustering(S_comb_all)

        S_comb_sel = S_comb_sel.assign(half = pd.Series(S_comb_all['half'],index=S_comb_sel.index).values)
        self.do_write_hdbscan_clustering_info(S_comb_sel)

    def do_write_hdbscan_clustering_info(self, S_comb_sel):

        ''' 
        Write clustering summary information 
        (i.e. cluster number, average scores, number of models)
        '''
        
        clusters = list(set(S_comb_sel['cluster']))
        
        out = open(self.analysis_dir+'summary_hdbscan_clustering.dat', 'w')
        out.write(' '.join(S_comb_sel.columns.values) +' n_models n_A n_B \n')
        for cl in clusters:
            cluster_info = np.array(S_comb_sel[S_comb_sel['cluster']==cl].mean().values).astype('str')
            n_models = len(S_comb_sel[S_comb_sel['cluster']==cl])
            n_A = len(S_comb_sel[(S_comb_sel['cluster']==cl) & (S_comb_sel['half']=='A')])
            n_B = len(S_comb_sel[(S_comb_sel['cluster']==cl) & (S_comb_sel['half']=='B')])
            out.write(' '.join(cluster_info.astype('str'))+' '+str(n_models)+' '+str(n_A)+' '+str(n_B) +'\n')
            print('Cluster number, number of models: ', cl, n_models)
        out.close()

    def do_write_hdbscan_clustering(self, S_comb):

        ''' 
        Write the frames information for each cluster
        '''

        print('Selecting and writing models to extract ...')
        clusters = list(set(S_comb['cluster']))
        clusters = [cl for cl in clusters if cl >=0]
        
        clus_sel = 0
        for cl in clusters:
            HH_cluster = S_comb[S_comb['cluster'] == cl]
            # Select two-halves
            HA = HH_cluster[(HH_cluster['half']=='A')]
            HB = HH_cluster[(HH_cluster['half']=='B')]
            if len(HA['half'].values) >= 10 and len(HB['half'].values) >= 10:
                clus_sel += 1
                n = self.plot_scores_distributions(HA, HB, cl)
            
                out_A = open(self.analysis_dir+'selected_models_A_cluster'+str(cl)+'.dat', 'w')
                out_B = open(self.analysis_dir+'selected_models_B_cluster'+str(cl)+'.dat', 'w')
                for row in HH_cluster.itertuples():
                    if row.half=='A':
                        out_A.write('h1_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3 \n')
                    if row.half=='B':
                        out_B.write('h2_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3 \n')

                out_A.close()
                out_B.close()

                out_A_det = open(self.analysis_dir+'selected_models_A_cluster'+str(cl)+'_detailed.dat', 'w')
                out_B_det = open(self.analysis_dir+'selected_models_B_cluster'+str(cl)+'_detailed.dat', 'w')
                for row in HH_cluster.itertuples():
                    if row.half=='A':
                        out_A_det.write('h1_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3 '+row.traj+' '+row.rmf3_file+' '+str(int(row.MC_frame))+' '+str(int(row.rmf_frame_index))+'\n')
                    if row.half=='B':
                        out_B_det.write('h2_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3 '+row.traj+' '+row.rmf3_file+' '+str(int(row.MC_frame))+' '+str(int(row.rmf_frame_index))+'\n')

                out_A_det.close()
                out_B_det.close()

                # Write to csv
                HA.to_csv(self.analysis_dir+'selected_models_A_cluster'+str(cl)+'_detailed.csv')
                HB.to_csv(self.analysis_dir+'selected_models_B_cluster'+str(cl)+'_detailed.csv')

                # Select n model from
                if int(n) < 30000 and len(HH_cluster) < 30000:
                    HH_sel = HH_cluster
                    HH_sel_A = HH_sel[(HH_sel['half']=='A')]
                    HH_sel_B = HH_sel[(HH_sel['half']=='B')]
                else:
                    HH_sel = HH_cluster.sample(n=29999)
                    HH_sel_A = HH_sel[(HH_sel['half']=='A')]
                    HH_sel_B = HH_sel[(HH_sel['half']=='B')]
                    
                HH_sel_A.to_csv(self.analysis_dir+'selected_models_A_cluster'+str(cl)+'_random.csv')
                HH_sel_B.to_csv(self.analysis_dir+'selected_models_B_cluster'+str(cl)+'_random.csv')
                
        if clus_sel == 0:
            print('WARNING: No models were selected because the simulations are not converged.')
                
    def plot_hdbscan_clustering(self, S_comb_sel, selected_scores):

        print('Generating HDBSCAN clustering plot ...')
        num_sel = len(S_comb_sel)-2

        palette = sns.color_palette("deep", len(np.unique(S_comb_sel['cluster']))).as_hex()
        cluster_colors = [palette[col] for col in S_comb_sel['cluster']]
    
        nn = len(selected_scores)-1
        fig, ax = pl.subplots(figsize=(5.0*nn, 5.0), nrows=1, ncols=nn)
        if nn == 1:
            ax = [ax]
        for i in range(len(selected_scores)-1):
            score_0 = selected_scores[i]
            score_1 = selected_scores[i+1]
            ax[i].scatter(S_comb_sel[score_0],S_comb_sel[score_1], c=np.array(cluster_colors),s=3.0,alpha=0.3)
            ax[i].set_title('HDBSCAN clustering')
            ax[i].set_xlabel(score_0 +' (a.u.)')
            ax[i].set_ylabel(score_1 +' (a.u.)')
            
        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)

        fig.savefig(self.analysis_dir+'plot_clustering_scores.png') 
        pl.close()
        
    def do_extract_models(self, gsms_info, filename, gsms_dir):

        # Split the DF
        df_array = np.array_split(gsms_info, self.nproc)

        # Define an output queue
        output = mp.Queue()
        
        # Setup a list of processes that we want to run
        processes = [mp.Process(target=self.extract_models, args=(df_array[x], filename, gsms_dir)) for x in range(self.nproc)]

        # Run processes
        for p in processes:
            p.start()
            
        # Exit the completed processes
        for p in processes:
            p.join() 
    
    def write_GSMs_info(self, gsms_info, filename):
        restraint_cols = gsms_info.columns.values
        gsms_info.to_csv(self.analysis_dir+filename, index=False)

    def get_models_to_extract(self, file):
        # Get models to extract from file
        DD  = pd.read_csv(file)
        return DD
        
    def get_sample_of_models_to_extract(self, file_A, file_B):
        DD_A = pd.read_csv(file_A)
        DD_B = pd.read_csv(file_B)

        # Take samples and see when it plateaus
        
    def extract_models(self, gsms_info, filename, gsms_dir):
        '''
        Use rmf_slice to extract the GSMs
        '''
        out_scores = open(gsms_dir+'.txt', 'w')
        for row in gsms_info.itertuples():
            id = row.traj
            fr = int(row.MC_frame)
            fr_rmf = int(row.rmf_frame_index)
            file = row.rmf3_file
            os.system('rmf_slice '+ id+'/'+file + ' '+gsms_dir+'/'+filename+'_'+str(id)+'_'+str(fr)+'.rmf3 --frame '+str(fr_rmf) )
            
            # Write score to file
            out_scores.write(str(row.Total_score)+'\n')
        out_scores.close()
            
    def create_gsms_dir(self, dir):
        '''
        Create directories for GSM.
        If already present, rename old one
        '''
        
        if not os.path.isdir(dir):
            os.makedirs(dir)
        else:
            os.system('mv '+dir + ' '+dir+'.old_'+str(random.randint(0,100)))
            os.makedirs(dir)

    def plot_scores_distributions(self, HA, HB, cl):
        '''
        Plot distribution of GSMs of both halves
        '''
        
        n_bins = 20
 
        scores_A = HA['Total_score']
        scores_B = HB['Total_score']

        fig = pl.figure(figsize=(8,4))
        gs = gridspec.GridSpec(1, 2,
                               width_ratios = [0.5,0.5],
                               height_ratios = [1.0])
    
        # Plot distributions
        ax= pl.subplot(gs[0])
        ax.hist(scores_A,  n_bins, histtype='step', stacked=True, fill=False, color='orangered')
        ax.hist(scores_B,  n_bins, histtype='step', stacked=True, fill=False, color='blue')
        ax.set_xlabel('Total scores (a.u.)')
        ax.set_ylabel('Density')

        # Plot expected score from random set
        nn = len(scores_A) + len(scores_B)
       
        M = np.arange(int(min(len(scores_A),len(scores_B))/20.), min(len(scores_A),len(scores_B)),int(nn/20.))
        if len(M)<=10.0:
            print(len(scores_A), len(scores_B))
            M = np.arange(int(min(len(scores_A),len(scores_B))/10.), min(len(scores_A),len(scores_B)),int(min(len(scores_A),len(scores_B))/10.))
    
        RH1 = []
        RH2 = []
        for m in M[1:]:
            D1 = []
            D2 = []
            for d in range(20):
                D1.append(min(random.sample(list(scores_A),m)))
                D2.append(min(random.sample(list(scores_B),m)))
            RH1.append((m,np.mean(D1),np.std(D1)))
            RH2.append((m,np.mean(D2),np.std(D2)))
        
        RH1 = np.array(RH1)
        RH2 = np.array(RH2)

        hits = 0
        for s in range(len(RH1)-1):
            dA = RH1[s,1]-RH1[s+1,1]
            dB = RH2[s,1]-RH2[s+1,1]
            if dA < RH1[s+1,2] and dB < RH2[s+1,1]:
                hits += 1
            if hits == 4:
                n = RH1[s+1,0]
                continue
        
        ax = pl.subplot(gs[1])
        ax.errorbar(RH1[:,0], RH1[:,1], yerr=RH1[:,2], c='orangered',fmt='o')
        ax.errorbar(RH2[:,0], RH2[:,1], yerr=RH2[:,2], c='blue',fmt='o')
        ax.axvline(n, color='grey')
        ax.set_xlim([0.8*M[0],1.1*M[-1]])
        ax.set_xlabel('Number of models')
        ax.set_ylabel('Minimum score (a.u.)')

        fig.savefig(self.analysis_dir+'plot_scores_convergence_cluster'+str(cl)+'.pdf') 
        pl.close()

        return n

    def get_XLs_details(self, XLs_type = None):
        '''
        For GSM, determine for each XLs how often it is satisfied.
        '''

        if len(self.S_dist_all) == 0:
            return 0
        
        traj_0 = self.S_dist_all.keys()[0]
        if XLs_type:
            col_sel = [v for v in self.S_dist_all[traj_0].columns.values if 'Distance' in v and XLs_type in v]
        else:
            col_sel = [v for v in self.S_dist_all[traj_0].columns.values if 'Distance' in v]
        i = 0
        # 1. Get distances for all models
        for traj in self.S_dist_all.keys():
            t = [x for x in traj.split('/') if self.dir_name in x][0]
            traj_number = int(t.split(self.dir_name)[1])
            T = self.S_dist_all[traj]
            T.to_csv(self.analysis_dir+'XLs_distances_'+str(traj_number)+'.csv')
            if i ==0 :
                all_dists = T
            else:
                all_dist.append(T, ignore_dist = True)
        
        # 2. All distances, from all traj
        all_dists.to_csv(self.analysis_dir+'XLs_distances_all.csv')
         
        # 3. Get XLs cutoff
        if XLs_type:
            cutoff = [v for k,v in self.XLs_cutoffs.items() if XLs_type in k][0]
        else:
            cutoff = list(self.XLs_cutoffs.values())[0]

        # 4. Check, per XLs, how often it is satisfied
        satif_XLs_indv =  all_dists.apply(lambda x: float(len(x[x<cutoff]))/float(len(x)), axis = 0)
        if XLs_type:
            satif_XLs_indv.to_csv(self.analysis_dir+'XLs_satisfaction_all_'+XLs_type+'.csv')    
        else:
            satif_XLs_indv.to_csv(self.analysis_dir+'XLs_satisfaction_all.csv')

       
    def plot_pEMAP_distances(self, pEMAP_satif, file_out):
        n_bins = 20
        fig, ax = pl.subplots(figsize=(4.0, 4.0), nrows=1, ncols=1)
        #axes[0].scatter(np.mean(S_dist[:,1:],axis=0), d_pEMAP, color='orangered',alpha=0.9)
        #xes[0].set_title('Dist pE-MAP', fontsize=14)
        #axes[0].set_xlabel('Dist. struct. (A)',fontsize=12)
        #axes[0].set_ylabel('Dist. model (A)',fontsize=12)
        
        ax.plot(pEMAP_satif[::10,0], pEMAP_satif[::10,1], color='orangered',alpha=0.8)
        ax.set_title('pE-MAP restraint', fontsize=14)
        ax.set_xlabel('Step',fontsize=12)
        ax.set_ylabel('Percent satisfied',fontsize=12)
        
        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(self.analysis_dir+file_out)

    def plot_Occams_satisfaction(self, Occams_satif, file_out):
        n_bins = 20
        fig, ax = pl.subplots(figsize=(4.0, 4.0), nrows=1, ncols=1)
        
        
        ax.plot(Occams_satif[::10,0], Occams_satif[::10,1], color='orangered',alpha=0.8)
        ax.set_title('Occams restraint', fontsize=14)
        ax.set_xlabel('Step',fontsize=12)
        ax.set_ylabel('Percent satisfied',fontsize=12)
        
        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(self.analysis_dir+file_out)
        
    def substrings(self, s):
        for i in range(len(s)):
            for j in range(i, len(s)):
                yield s[i:j+1]
        
    def get_str_match(self, strs):
        if len(strs) > 1:
            intersect = set(self.substrings(strs[0])) & set(self.substrings(strs[1]))
            return max(intersect, key = len)
        else:
            return strs[0]


        
