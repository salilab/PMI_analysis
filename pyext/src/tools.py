#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   Tools for analysis of MCMC
"""
from __future__ import print_function
import os
import glob


class ReadClustering(object):
    '''Read output from structural clustering'''
    def __init__(self,
                 clustering_dir):

        self.clustering_dir = clustering_dir

    def get_number_of_clusters(self):

        clusters = glob.glob(os.path.join(
            self.clustering_dir, 'cluster.*.all.txt'))
        return len(clusters)

    def get_rmfs_cluster(self, cluster):
        '''
        Read identities of model in clusters and then
        get the names of all rmfs
        '''

        rmfs_dic = {}
        print('------', self.clustering_dir)
        with open(os.path.join(self.clustering_dir, 'Identities_A.txt'), 'r') as f:
            for line in f:
                vals = line.split()
                rmfs_dic[vals[1]] = vals[0]
        with open(os.path.join(self.clustering_dir, 'Identities_B.txt'), 'r') as f:
            for line in f:
                vals = line.split()
                rmfs_dic[vals[1]] = vals[0]

        # Read rmfs in cluster
        rmfs = []
        with open(os.path.join(
                self.clustering_dir,
                'cluster.'+str(cluster)+'.sample_A.txt'), 'r') as f:
            for mod in f:
                vals = mod.split()[0]
                if vals in rmfs_dic:
                    rmfs.append(rmfs_dic[vals])
                else:
                    print('Model missing: ', vals)
        with open(os.path.join(
                self.clustering_dir,
                'cluster.'+str(cluster)+'.sample_B.txt'), 'r') as f:
            for mod in f:
                vals = mod.split()[0]
                if vals in rmfs_dic:
                    rmfs.append(rmfs_dic[vals])
                else:
                    print('Model missing: ', vals)

        return rmfs
