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

import IMP
import IMP.rmf
import RMF

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})


class ValidationModels(object):
    def __init__(self,
                 clustering_dir,
                 scores_sample_A,
                 scores_sample_B):

        self.clustering_dir = clustering_dir
        self.scores_sample_A = scores_sample_A
        self.scores_sample_B = scores_sample_B
        
        self.manager = mp.Manager()
        
        self.read_scores_files()
        self.read_identities()
        self.read_clusters()
        self.read_clustered_models()
        
    def read_scores_files(self):
        # Read scores from file
        S1 = pd.read_csv(self.scores_sample_A, sep=',')
        S2 = pd.read_csv(self.scores_sample_B ,sep=',')
        self.S = pd.concat([S1,S2])

    def read_identities(self):
        
        # Read identities of models in each cluster
        ids_file_A = pd.read_csv(self.clustering_dir+'/Identities_A.txt', sep=' ', names=['model','id'])
        ids_file_B = pd.read_csv(self.clustering_dir+'/Identities_B.txt', sep=' ', names=['model','id'])

        ids_file_A['file_name'] = ids_file_A['rmf3_file'].apply(lambda model: model.split('/')[-1])
        ids_file_B['file_name'] = ids_file_B['rmf3_file'].apply(lambda model: model.split('/')[-1])
        ids_file_A['half'] = 'A'
        ids_file_B['half'] = 'B'
        
        self.DC = pd.concat([ids_file_A, ids_file_B])
        self.DC['traj'] = np.nan
        self.DC['MC_frame'] = np.nan
        self.DC['cluster'] = np.nan

        print(self.DC.head())
        
    def read_clusters(self):
        '''
        Get all rmf3 files in clusters
        '''
    
        self.cluster_rmfs = {}

        clusters = glob.glob(self.clustering_dir+'/cluster.*.all.txt').sort()
        self.n_cluster = len(clusters)
        for n, cl in enumerate(clusters):
            print(n, cl)
            with open(cluster) as fh:
                ids = f.readlines()
            print(ids)
            self.DC[self.DC[id].isin(ids)]['cluster'] = n
            
        print(self.DC.head())
            
    def read_clustered_models(self):
        '''
        For each cluster, get the models information
        and get worst scoring model for restraints
        '''
        
        score_fields = [s for s in self.S.columns.values if 'sum' in s]
        self.sel_frames = {}
        
        for n range(self.n_clusters):
            max_scores = {s: [0.0,''] for s in score_fields}

            frames = self.DC[self.DC['cluster']==n]['file_name']
            print(frames.head())
            
            idx = pd.Index(self.S['frame_RMF3']).get_indexer(frames.values)

            info_scores = self.S.iloc[idx]
                
            self.DC[self.DC['cluster']==n]['traj'] = info_scores['traj'].values
            self.DC[self.DC['cluster']==n]['MC_frame'] = info_scores['MC_frame'].values
            
            print(self.DC.head())

            for s in score_fields:
                m = info_scores[s].max()
                f = info_scores[info_scores[s]==m]['frame_RMF3'].values[0]
                max_scores[s] = [m,f]
                
            self.sel_frames[n] = max_scores

        print(self.sel_frames)
                
    def excluded_volume_satisfaction(self,
                                     rmf3_file,
                                     resolution=10):
        '''
        Take the worst scoring model in the cluster and 
        compute the EV satisfaction 
        '''

        kappa = 1.0
        
    
        m = IMP.Model()
        hier = IMP.pmi.analysis.get_hiers_from_rmf(m,0,rmf3_file)
        
        included_ps = IMP.atom.Selection(hier).get_selected_particles()
        ppis = {p.get_index() : p for p in included_ps}
           
        lsa = IMP.container.ListSingletonContainer(m)
        lsa.add(IMP.get_indexes(included_ps))

        rbcpf = IMP.core.RigidClosePairsFinder()
        cpc = IMP.container.ClosePairContainer(lsa, 0.0, rbcpf, 10.0)
        
        cps = rbcpf.get_close_pairs(m, included_ps)
        D = []
        for i0, i1 in cps:            
            c0 = IMP.core.XYZR(ppis[i0])
            c1 = IMP.core.XYZR(ppis[i1])
            dist=IMP.core.get_distance(c0, c1)
            D.append(dist)

        viol = [v for v in D if v < -kappa]

        print(100.*(1-len(viol)/float(len(D))))

        return round(100.*(1-len(viol)/float(len(D))),2)
        
    def sequence_connectivity_satisfaction(self, rmf3_file):

        # TO DO!

        return 0

    def get_XLs_distances(self):

       self.XLs_dist_clusters = {} 
        
        for n in range(self.n_clusters):
            dist_all = pd.DataFrame()
            sel_cluster = self.DC[self.DC['cluster']==n]
            trajs =  sel_cluster['traj'].apply(lambda x: x.split('run_')[1]).unique()
            for t in traj:
                frames = self_cluster['MC_frame']
                dist = pd.read_csv(self.analys_dir+'/XLs_dist_info_'+str(t)+'.csv')
                if not dist_all.empty:
                    dist_all.append(dist)
                else:
                    dist_all = dist

            self.XLs_dist_clusters[n] = dist_all
                    

    def summarize_XLs_info(self):

        if self.Multiple_XLs_restraints:
            columns = ['cluster'] + \
                      [t +'_ensemble_satisfaction' for t in self.XLs_cutoffs.keys()] + \
                      [t +'_average_frame_satisfaction' for t in self.XLs_cutoffs.keys()]
        else:
            columns = ['cluster','ensemble_satisfaction','average_frame_satisfaction']
        
        sXLs = pd.DataFrame(names = columns)
        
        for n in range(self.n_clusters):
            l = [n]
            if self.Multiple_XLs_restraints:
                for type_XLs in self.XLs_cutoffs.keys():
                    cutoff = self.XLs_cutoffs[type_XLs]
                    ensemble, ave = get_XLs_statistics(self,
                                                       cluster=n,
                                                       type_XLs = type_XLs,
                                                       cutoff=cutoff)
                    l = l + [ensemble+ave]
                sXLs.iloc[n] = l
                    
            else:
                cutoff = list(self.XLs_cutoffs.values())[0]
                ensemble, ave = self.get_XLs_details(cluster = cl,
                                                      cutoff=cutoff)
                sXLs.iloc[n] = [n, ensemble, ave]

        sXLs.to_csv(self.clustering_dir+'XLs_validation_summary.csv')
            
        
    def get_XLs_statistics(self, cluster=0, type_XLs = None, cutoff=30.):
        '''
        For each cluster, determine for each XLs, how often it is satisfied.
        '''
        
        XLs_cluster =  self.XLs_dist_clusters[n]
        
        if type_XLs:
            dist_columns = [v for v in XLs_cluster.columns.values if 'Distance' in v and type_XLs in v]
            cutoff = [v for k,v in self.XLs_cutoffs.items() if type_XLs in k][0]
        else:
            dist_columns = [v for v in XLs_cluster.columns.values if 'Distance' in v]
            cutoff = list(self.XLs_cutoffs.values())[0]

    
        dXLs_cluster = XLs_cluster.loc[:, dist_columns]
            
        stats_XLs = pd.DataFrame()
        stats_XLs['mean'] = dXLs_cluster.mean()
        stats_XLs['std'] = dXLs_cluster.std()
        stats_XLs['min'] = dXLs_cluster.min()
        stats_XLs['max'] = dXLs_cluster.max()
        
        stats_XLs['perc_satif'] = dXLs_cluster.apply(lambda x: float(len(x[x<cutoff]))/float(len(x)), axis = 0)

        if type_XLs:
            stats_XLs.to_csv(self.clustering_dir+'XLs_satisfaction_'+type_XLs+'_cluster_'+str(cluster)+'.csv')
            dXLs_cluster.to_csv(self.clustering_dir+'XLs_distances_'+type_XLs+'_cluster_'+str(cluster)+'.csv')
        else:
            stats_XLs.to_csv(self.clustering_dir+'XLs_satisfaction_cluster_'+str(cluster)+'.csv')
            dXLs_cluster.to_csv(self.clustering_dir+'XLs_distances_cluster_'+str(cluster)+'.csv')

        # Now compute the frame and ensemble satisfaction
        ensemble = dXLs_cluster.apply(lambda x: float(len(x[x>cutoff]))/len(x)) ,axis=0).sum()
        ave = dXLs_cluster.apply(lambda x: float(len(x[x<cutoff]))/len(x) ,axis=1).mean()
        
        return ensemble, ave
            
    def XLs_satisfaction_from_dataset(self, XLs_dataset = None):

        # TO DO

        return 0


    def get_clusters_info(self):

        self.info_clusters = {} 
        
        for n in range(self.n_clusters):
            info_all = pd.DataFrame()
            sel_cluster = self.DC[self.DC['cluster']==n]
            trajs =  sel_cluster['traj'].apply(lambda x: x.split('run_')[1]).unique()
            for t in traj:
                frames = self_cluster['MC_frame']
                info = pd.read_csv(self.analys_dir+'/other_info_'+str(t)+'.csv')
                if not info_all.empty:
                    info_all.append(info)
                else:
                    info_all = info

            self.info_clusters[n] = info_all
    
    def pEMAP_satisfaction(self):

        sPEMAP = pd.DataFrame(names = ['cluster','pEMAP_satisfaction', 'pEMAP_sigma'])
        
        for n in range(self.n_clusters):
            m = self.info_clusters[n]['pEMapRestraint_satisfaction'].mean()
            sigma = self.info_clusters[n]['pEMapRestraint_sigma_0'].mean()

            sPEMAP.iloc[n] = [n, m, sigma]

        sPEMAP.to_csv(self.clustering_dir+'pEMAP_validation_summary.csv')

    def EM3D_satisfaction(self):

        # TO DO

        return 0



