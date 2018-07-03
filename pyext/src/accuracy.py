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

class AccuracyModels(object):
    def __init__(self,
                 selection_dictionary,
                 clustering_dir,
                 scores_sample_A,
                 scores_sample_B,
                 ref_rmf3,
                 nproc=6):

        self.selection_dictionary = selection_dictionary
        self.clustering_dir = clustering_dir
        self.scores_sample_A = scores_sample_A
        self.scores_sample_B = scores_sample_B
        self.ref_rmf3 = ref_rmf3       

        self.manager = mp.Manager()

        self.all_accu = {}
        self.ids_all = {}
 
        # Read all models files
        self.read_scores_files()

        # Read cluster files
        clusters = glob.glob(self.clustering_dir+'cluster.*.all.txt')
        n_clusters = len(clusters)

        self.all_accu = {}
        for cl in clusters:
            cl_number = cl.split('.')[-3]
            rmfs_clus = []
            for line in open(cl, 'r'):
                rmf_file = self.ids_all[line.split('\n')[0]]
                rmfs_clus.append(rmf_file)
            
            rmfs_clus = rmfs_clus[::2]
            # Divide rmfs into groups
            ND = int(np.ceil(len(rmfs_clus)/float(nproc)))
            rmfs_dict = {}
            for k in range(nproc-1):
                rmfs_dict[k] = rmfs_clus[(k*ND):(k*ND+ND)]
            rmfs_dict[nproc-1] = rmfs_clus[((nproc-1)*ND):(len(rmfs_clus))]

            # Compute accuracies
            cl_accu = self.manager.list()
            # Define an output queue
            output = mp.Queue()
        
            # Setup a list of processes that we want to run
            processes = [mp.Process(target=self.accuracy_rmfs, args=(rmfs_dict[x], cl_accu)) for x in range(nproc)]

            # Run processes
            for p in processes:
                p.start()

            # Exit the completed processes
            for p in processes:
                p.join()

            # Add to dictionary
            self.all_accu[cl_number] = cl_accu

        # Plot and write results
        self.plot_accuracy_histograms()
        self.write_accuracy_values()
        self.plot_score_versus_accuracy()        

    def accuracy_rmfs(self, rmfs, cl_accu):
        if len(rmfs)>0:
            frames=[0]*len(rmfs)
            model=IMP.Model()
            pr=IMP.pmi.analysis.Precision(model,
                                          resolution=1,
                                          selection_dictionary=self.selection_dictionary)
            pr.set_precision_style('pairwise_rmsd')
            pr.add_structures(zip(rmfs,frames), 'set0')
            pr.set_reference_structure(self.ref_rmf3,0)
            vals = pr.get_rmsd_wrt_reference_structure_with_alignment('set0',('all_sel'))
            cl_accu.extend(zip(rmfs, vals['all_sel']['all_distances']))

    def read_scores_files(self):
        # Read scores from file
        S1 = pd.read_table(self.scores_sample_A, sep=',')
        S2 = pd.read_table(self.scores_sample_B ,sep=',')
        self.S = pd.concat([S1,S2])
        
        # Read identities of models, and put into dictionary
        ids_B = {}
        ids_file_A = open(self.clustering_dir+'/Identities_A.txt','r')
        ids_file_B = open(self.clustering_dir+'/Identities_B.txt','r')
    
        for line in ids_file_A:
            vals = line.split()
            self.ids_all[vals[1]] = vals[0]
        for line in ids_file_B:
            vals = line.split()
            self.ids_all[vals[1]] = vals[0]

    def write_accuracy_values(self):
        self.score_accu = {}
        out_summary = open(self.clustering_dir+'/accuracy_summary_clusters.dat','w')
        out_summary.write('cluster mean min max N \n')
        for k, v in self.all_accu.iteritems():
            out = open(self.clustering_dir+'/accuracy_cl'+str(k)+'.dat','w')
            vv = np.array(v)[:,1].astype(float)
            out_summary.write(str(k)+'\t'+str(np.mean(vv))+'\t'+str(np.min(vv))+'\t'+str(np.max(vv))+'\t'+str(len(vv))+'\n')
            T = []
            for f, accu in v:
                f_short = f.split('/')[-1]
                run = 'run_'+f_short.split('_')[2]
                frame = f_short.split('_')[3].split('.rmf3')[0]
                score = self.S[(self.S['traj'] == run) & (self.S['MC_frame'] == float(frame))]['Total_score'].values[0]
                T.append([score, accu])
                out.write('%s\t%s\t%s \n'%(f_short,accu,score))
            out.close()
            self.score_accu[k] = np.array(T)
        out_summary.close()

    def plot_accuracy_histograms(self):
        n_bins = 10
        clus_all = sorted(self.all_accu.keys())
        fig, ax = pl.subplots(figsize=(5.0, 5.0), nrows=1, ncols=1)
        palette = pl.get_cmap('Set1')
        colors = [palette(1.*i/len(clus_all)) for i in range(len(clus_all))]
        for i, clus in enumerate(clus_all):
            A = np.array(self.all_accu[clus])[:,1].astype(float)
            ax.hist(A, n_bins,  histtype='step',fill=False, color=colors[i],alpha=0.9)
            ax.axvline(np.mean(A), color=colors[i], alpha=0.9)
             
        ax.set_title('Accuracy of selected models', fontsize=14)
        ax.set_xlabel('Accuracy (A)',fontsize=12)
        ax.set_ylabel('Number of models',fontsize=12)

        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(self.clustering_dir+'/plot_accuracy_all_clusters.pdf')
        
        
    def plot_score_versus_accuracy(self):
        # Plot score vs accuracy for first cluster
        fig, ax = pl.subplots(figsize=(5.0, 5.0), nrows=1, ncols=1)
        ax.scatter(self.score_accu['0'][:,0],self.score_accu['0'][:,1],color='royalblue', alpha=0.8, s=30, edgecolors='none')
        ax.set_title('Score vs. Accuracy of selected models (cluster 0)', fontsize=12)
        ax.set_xlabel('Total Score (a.u.)',fontsize=12)
        ax.set_ylabel('Accuracy (A)',fontsize=12)
        
        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(self.clustering_dir+'/plot_score_vs_accuracy.pdf')
