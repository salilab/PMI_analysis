#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   Tools for analysis of MCMC
"""
from __future__ import print_function
import IMP
import IMP.algebra
import IMP.em
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.output
import IMP.rmf
import RMF
import IMP.pmi.analysis
from operator import itemgetter
from copy import deepcopy
from math import log,sqrt
import itertools
import numpy as np
import scipy.spatial.distance
import os
import math
import glob

import multiprocessing as mp

class ReadClustering(object):
    '''Read output from structural clustering'''
    def __init__(self,
                 clustering_dir):

        self.clustering_dir = clustering_dir

    def get_number_of_clusters(self):

        clusters = glob.glob(os.path.join(
            self.clustering_dir,'cluster.*.all.txt'))
        return len(clusters)
        
    def get_rmfs_cluster(self, cluster):
        '''
        Read identities of model in clusters and then
        get the names of all rmfs
        '''

        rmfs_dic = {}
        print('------', self.clustering_dir)
        for line in open(os.path.join(self.clustering_dir,'Identities_A.txt'), 'r'):
            vals = line.split()
            rmfs_dic[vals[1]] = vals[0]
        for line in open(os.path.join(self.clustering_dir,'Identities_B.txt'), 'r'):
            vals = line.split()
            rmfs_dic[vals[1]] = vals[0]

        # Read rmfs in cluster
        rmfs = []
        for mod in open(os.path.join(
                self.clustering_dir,'cluster.'+str(cluster)+'.sample_A.txt'),'r'):
            vals = mod.split()[0]
            try:
                rmfs.append(rmfs_dic[vals])
            except:
                print('Model missing: ', vals)
        for mod in open(os.path.join(
                self.clustering_dir,'cluster.'+str(cluster)+'.sample_B.txt'),'r'):
            vals = mod.split()[0]
            try:
                rmfs.append(rmfs_dic[vals])
            except:
                print('Model missing: ', vals)

        return rmfs
