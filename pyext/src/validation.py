from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.pmi.tools
import IMP.atom
import glob
import math
import numpy as np
import pandas as pd
import itertools
import os
import multiprocessing as mp

import IMP
import IMP.rmf

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 8})

import tools  # noqa: E402


class ValidationModels(object):
    def __init__(self,
                 analysis_dir,
                 clustering_dir,
                 scores_sample_A,
                 scores_sample_B,
                 XLs_cutoffs):

        self.analysis_dir = analysis_dir
        self.clustering_dir = clustering_dir
        self.scores_sample_A = scores_sample_A
        self.scores_sample_B = scores_sample_B
        self.XLs_cutoffs = XLs_cutoffs

        print(self.clustering_dir)

        self.manager = mp.Manager()

        self.read_scores_files()
        self.read_identities()
        self.read_clusters()
        self.read_clustered_models()

    def read_scores_files(self):
        # Read scores from file
        S1 = pd.read_csv(self.scores_sample_A, sep=',')
        S2 = pd.read_csv(self.scores_sample_B, sep=',')
        self.S = pd.concat([S1, S2])

        if 'frame_RMF3' not in self.S.columns:
            self.S.loc[:, 'frame_RMF3'] = self.S.apply(
                lambda row: 'h1_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3'
                if row.half == 'A'
                else 'h2_'+row.traj+'_'+str(int(row.MC_frame))+'.rmf3',
                axis=1)

    def read_identities(self):

        # Read identities of models in each cluster
        ids_file_A = pd.read_csv(
            os.path.join(self.clustering_dir, 'Identities_A.txt'),
            sep=' ', names=['model', 'id'])
        ids_file_B = pd.read_csv(
            os.path.join(self.clustering_dir, 'Identities_B.txt'),
            sep=' ', names=['model', 'id'])

        ids_file_A['file_name'] \
            = ids_file_A['model'].apply(lambda model: model.split('/')[-1])
        ids_file_B['file_name'] \
            = ids_file_B['model'].apply(lambda model: model.split('/')[-1])
        ids_file_A['half'] = 'A'
        ids_file_B['half'] = 'B'

        self.DC = pd.concat([ids_file_A, ids_file_B])
        self.DC['traj'] = np.nan
        self.DC['MC_frame'] = np.nan
        self.DC['cluster'] = np.nan

    def read_clusters(self):
        '''
        Get all rmf3 files in clusters
        '''

        self.cluster_rmfs = {}

        clusters = np.sort(
            glob.glob(os.path.join(self.clustering_dir, 'cluster.*.all.txt')))
        self.n_clusters = len(clusters)
        for n, cluster in enumerate(clusters):
            ids = pd.read_csv(cluster, names=['id'])
            v = ids['id'].values
            self.DC.loc[self.DC['id'].isin(v), 'cluster'] = int(n)

    def read_clustered_models(self):
        '''
        For each cluster, get the models information
        and get worst scoring model for restraints
        '''

        score_fields = [s for s in self.S.columns.values if 'sum' in s]
        self.sel_frames = {}

        for n in range(self.n_clusters):
            max_scores = {s: [0.0, ''] for s in score_fields}

            frames = self.DC[self.DC['cluster'] == n]['file_name']
            print(self.S.head())
            for a in self.S.columns:
                print(a)
            print(self.S['rmf3_file'].head())
            idx = pd.Index(self.S['frame_RMF3']).get_indexer(frames.values)

            info_scores = self.S.iloc[idx]

            self.DC.loc[self.DC['cluster'] == n, 'traj'] \
                = info_scores['traj'].values
            self.DC.loc[self.DC['cluster'] == n, 'MC_frame'] \
                = info_scores['MC_frame'].values

            for s in score_fields:
                m = info_scores[s].max()
                f = info_scores[info_scores[s] == m]['frame_RMF3'].values[0]
                max_scores[s] = [m, f]

            self.sel_frames[n] = max_scores

    def get_excluded_volume_satisfaction(self):

        sEV = pd.DataFrame(columns=['cluster', 'EV_satisfaction'])

        for n in range(self.n_clusters):
            rmf3_file = self.sel_frames[n]['EV_sum'][1]
            rmf3_full = self.DC[
                self.DC['file_name'] == rmf3_file]['model'].values[0]
            s = self.excluded_volume_satisfaction(rmf3_full)
            sEV.loc[n] = [n, s]

        sEV.to_csv(os.path.join(self.clustering_dir,
                                'EV_validation_summary.csv'),
                   index=False)

    def excluded_volume_satisfaction(self,
                                     rmf3_file,
                                     resolution=10):
        '''
        Take the worst scoring model in the cluster and
        compute the EV satisfaction
        '''

        kappa = 1.0

        m = IMP.Model()
        hier = IMP.pmi.analysis.get_hiers_from_rmf(m, 0, rmf3_file)

        included_ps = IMP.atom.Selection(hier).get_selected_particles()
        ppis = {p.get_index(): p for p in included_ps}

        lsa = IMP.container.ListSingletonContainer(m)
        lsa.add(IMP.get_indexes(included_ps))

        rbcpf = IMP.core.RigidClosePairsFinder()
        # cpc = IMP.container.ClosePairContainer(lsa, 0.0, rbcpf, 10.0)

        cps = rbcpf.get_close_pairs(m, included_ps)
        D = []
        for i0, i1 in cps:
            c0 = IMP.core.XYZR(ppis[i0])
            c1 = IMP.core.XYZR(ppis[i1])
            dist = IMP.core.get_distance(c0, c1)
            D.append(dist)

        viol = [v for v in D if v < -kappa]

        return round(100.*(1-len(viol)/float(len(D))), 2)

    def sequence_connectivity_satisfaction(self, rmf3_file):

        # TO DO!

        return 0

    def get_XLs_distances(self):

        self.XLs_dist_clusters = {}

        for n in range(self.n_clusters):
            dist_all = pd.DataFrame()
            sel_cluster = self.DC[self.DC['cluster'] == n]
            trajs = sel_cluster['traj'].apply(
                lambda x: x.split('run_')[1]).unique()
            for t in trajs:
                frames = sel_cluster[
                    sel_cluster['traj'] == 'run_'+t]['MC_frame']
                dist = pd.read_csv(self.analysis_dir + '/XLs_dist_info_'
                                   + str(t) + '.csv')
                dist_cluster = dist[dist['MC_frame'].isin(frames)]
                if not dist_all.empty:
                    dist_all.append(dist_cluster)
                else:
                    dist_all = dist_cluster

            self.XLs_dist_clusters[n] = dist_all

    def get_XLs_satisfaction(self):

        self.get_XLs_distances()
        self.Multiple_XLs_restraints = None

        if self.Multiple_XLs_restraints:
            columns = ['cluster'] + \
                      [t + '_ensemble_satisfaction'
                       for t in self.XLs_cutoffs.keys()] + \
                      [t + '_average_frame_satisfaction'
                       for t in self.XLs_cutoffs.keys()]
        else:
            columns = ['cluster', 'ensemble_satisfaction',
                       'average_frame_satisfaction']

        sXLs = pd.DataFrame(columns=columns)

        for n in range(self.n_clusters):
            iloc = [n]
            if self.Multiple_XLs_restraints:
                for type_XLs in self.XLs_cutoffs.keys():
                    cutoff = self.XLs_cutoffs[type_XLs]
                    ensemble, ave = self.XLs_statistics(cluster=n,
                                                        type_XLs=type_XLs,
                                                        cutoff=cutoff)
                    iloc += [ensemble+ave]
                sXLs.iloc[n] = iloc

            else:
                cutoff = list(self.XLs_cutoffs.values())[0]
                ensemble, ave = self.XLs_statistics(cluster=n, cutoff=cutoff)
                sXLs.loc[n] = [n, ensemble, ave]

        sXLs.to_csv(os.path.join(self.clustering_dir,
                                 'XLs_validation_summary.csv'),
                    index=False)

    def XLs_statistics(self, cluster=0, type_XLs=None, cutoff=30.):
        '''
        For each cluster, determine for each XLs, how often it is satisfied.
        '''

        XLs_cluster = self.XLs_dist_clusters[cluster]

        if type_XLs:
            dist_columns = [v for v in XLs_cluster.columns.values
                            if 'Distance' in v and type_XLs in v]
            cutoff = [v for k, v in self.XLs_cutoffs.items()
                      if type_XLs in k][0]
        else:
            dist_columns = [v for v in XLs_cluster.columns.values
                            if 'Distance' in v]
            cutoff = list(self.XLs_cutoffs.values())[0]

        dXLs_cluster = XLs_cluster.loc[:, dist_columns]

        stats_XLs = pd.DataFrame()
        stats_XLs['mean'] = dXLs_cluster.mean()
        stats_XLs['std'] = dXLs_cluster.std()
        stats_XLs['min'] = dXLs_cluster.min()
        stats_XLs['max'] = dXLs_cluster.max()

        stats_XLs['perc_satif'] = dXLs_cluster.apply(
            lambda x: float(len(x[x < cutoff]))/float(len(x)), axis=0)

        if type_XLs:
            stats_XLs.to_csv(
                os.path.join(self.clustering_dir,
                             'XLs_satisfaction_' + type_XLs + '_cluster_'
                             + str(cluster) + '.csv'))
            dXLs_cluster.to_csv(
                os.path.join(self.clustering_dir,
                             'XLs_distances_' + type_XLs + '_cluster_'
                             + str(cluster) + '.csv'))
        else:
            stats_XLs.to_csv(
                os.path.join(self.clustering_dir,
                             '/XLs_satisfaction_cluster_'+str(cluster)+'.csv'))
            dXLs_cluster.to_csv(
                os.path.join(self.clustering_dir,
                             '/XLs_distances_cluster_'+str(cluster)+'.csv'))

        # Now compute the frame and ensemble satisfaction
        min_ensemble = dXLs_cluster.apply(lambda x: min(x), axis=0)
        ensemble = len(min_ensemble[min_ensemble < cutoff]) \
            / float(len(min_ensemble))
        ave = dXLs_cluster.apply(
            lambda x: float(len(x[x < cutoff]))/len(x), axis=1).mean()

        return ensemble, ave

    def XLs_satisfaction_from_dataset(self, XLs_dataset=None):

        # TO DO

        return 0

    def get_clusters_info(self):

        self.info_clusters = {}

        for n in range(self.n_clusters):
            info_all = pd.DataFrame()
            sel_cluster = self.DC[self.DC['cluster'] == n]
            trajs = sel_cluster['traj'].apply(
                lambda x: x.split('run_')[1]).unique()
            for t in trajs:
                frames = \
                    sel_cluster[sel_cluster['traj'] == 'run_'+t]['MC_frame']
                info = pd.read_csv(os.path.join(
                    self.analysis_dir, 'other_info_'+str(t)+'.csv'))
                info_cluster = info[info['MC_frame'].isin(frames)]
                if not info_all.empty:
                    info_all.append(info_cluster)
                else:
                    info_all = info_cluster

            self.info_clusters[n] = info_all

    def get_pEMAP_satisfaction(self):

        self.get_clusters_info()

        sPEMAP = pd.DataFrame(columns=['cluster', 'pEMAP_satisfaction',
                                       'pEMAP_sigma'])

        for n in range(self.n_clusters):
            m = self.info_clusters[n]['pEMapRestraint_satisfaction'].mean()
            sigma = self.info_clusters[n]['pEMapRestraint_sigma_0'].mean()

            sPEMAP.loc[n] = [n, m, sigma]

        sPEMAP.to_csv(os.path.join(
            self.clustering_dir, 'pEMAP_validation_summary.csv'), index=False)

    def distance_implied_by_MIC(self, MIC_file, rmf):
        k = -0.014744
        n = -0.41
        dists_MIC = []

        model = IMP.Model()
        hier = IMP.pmi.analysis.get_hiers_from_rmf(model, 0, rmf)[0]

        # Check if residue pair is in representaiton
        for i, row in self.MIC.iterrows():
            s1 = IMP.atom.Selection(
                hier, molecule=row['p1'],
                residue_index=int(row['r1'])).get_selected_particles()
            s2 = IMP.atom.Selection(
                hier, molecule=row['p2'],
                residue_index=int(row['r2'])).get_selected_particles()

            if len(s1) > 0 and len(s2) > 0:
                m = row['MIC']
                if m <= 0.6:
                    d_mic = (math.log(0.6)-n)/k
                else:
                    d_mic = (math.log(m)-n)/k
                dists_MIC.append(d_mic)

        return dists_MIC

    def get_pEMAP_satisfaction_full(self, MIC_file):
        self.nproc = 10
        sPEMAP_all = pd.DataFrame(columns=['cluster', 'pEMAP_satisfaction'])

        # Read MIC file
        self.MIC = pd.read_csv(MIC_file, sep=' ',
                               names=['p1', 'p2', 'r1', 'r2', 'MIC', 'd'])
        print(self.MIC.head())

        # For each cluster read all rmfs and accumulated distances
        RC = tools.ReadClustering(self.clustering_dir)
        n_clusters = RC.get_number_of_clusters()

        for cl in range(n_clusters):
            self.dists_pEMAP = self.manager.dict()
            rmfs_cluster = RC.get_rmfs_cluster(cl)

            # Divide rmfs into groups
            ND = int(np.ceil(len(rmfs_cluster)/float(self.nproc)))
            rmfs_dict = {}
            for k in range(self.nproc-1):
                rmfs_dict[k] = rmfs_cluster[(k*ND):(k*ND+ND)]
            rmfs_dict[self.nproc-1] \
                = rmfs_cluster[((self.nproc-1)*ND):(len(rmfs_cluster))]

            # Setup a list of processes that we want to run
            processes = [mp.Process(target=self.distances_pEMAP,
                                    args=(rmfs_dict[x], 0))
                         for x in range(self.nproc)]

            # Run processes
            for p in processes:
                p.start()

            # Exit the completed processes
            for p in processes:
                p.join()

            # Save validation report
            dist_mic = self.distance_implied_by_MIC(MIC_file, rmfs_cluster[0])
            AA = [v for k, v in self.dists_pEMAP.items()]
            AA = np.array(AA)
            dist_min = np.apply_along_axis(np.min, axis=0, arr=AA)
            satif = [1 for i, j in zip(dist_min, dist_mic) if i >= j]
            percent_satif = float(len(satif))/len(dist_mic)
            sPEMAP_all = sPEMAP_all.append(
                {'cluster': cl, 'pEMAP_satisfaction': percent_satif},
                ignore_index=True)

        print(n_clusters, sPEMAP_all)
        sPEMAP_all.to_csv(
            os.path.join(self.clustering_dir,
                         'pEMAP_all_validation_summary.csv'),
            index=False)

    def distances_pEMAP(self, rmfs, t):

        for k, rmf in enumerate(rmfs):
            model = IMP.Model()
            hier = IMP.pmi.analysis.get_hiers_from_rmf(model, 0, rmf)[0]
            dists = self.get_all_distances_pEMAP(hier)
            self.dists_pEMAP[rmf] = dists
            del model, hier

    def get_distance_pair(self, s1, s2):
        dd = []
        for ss1, ss2 in itertools.product(s1, s2):
            c1 = IMP.core.XYZ(ss1)
            c2 = IMP.core.XYZ(ss2)
            dist = IMP.core.get_distance(c1, c2)
            dd.append(dist)
        return np.min(dd)

    def get_all_distances_pEMAP(self, hier):
        dist = []
        for i, row in self.MIC.iterrows():
            s1 = IMP.atom.Selection(
                hier, molecule=row['p1'],
                residue_index=int(row['r1'])).get_selected_particles()
            s2 = IMP.atom.Selection(
                hier, molecule=row['p2'],
                residue_index=int(row['r2'])).get_selected_particles()

            if len(s1) > 0 and len(s2) > 0:
                d = self.get_distance_pair(s1, s2)
                dist.append(d)

        return dist

    def EM3D_satisfaction(self):

        # TO DO

        return 0
