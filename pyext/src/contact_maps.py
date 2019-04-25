#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   Tools for contact maps analysis
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

import itertools
import numpy as np
import scipy.spatial.distance
import os
import math
import random

import multiprocessing as mp
import tools

class CMTable(object):
        ''' Compute contacts maps for all models in ensemble '''
        def __init__(self,
                     GSMs_dir = '',
                     clustering_dir = '',
                     out_dir = 'CMs/',
                     cluster = 0,
                     number_of_models = 10000,
                     selection = [],
                     cutoff = 16.0,
                     XLs_cutoff = 35.0,
                     nproc = 10):

                self.GSMs_dir = GSMs_dir
                self.clustering_dir = clustering_dir
                self.cluster = cluster
                self.number_of_models = number_of_models
                self.selection = selection
                self.cutoff = cutoff
                self.XLs_cutoff = XLs_cutoff
                self.nproc = nproc

                # Contact maps directory
                self.out_dir = out_dir
                if os.path.exists(self.out_dir):
                        print('Overwriting '+ self.out_dir)
                else:
                        os.makedirs(self.out_dir)
        
                # Dictionaries for results
                self.manager = mp.Manager()
                self.cm_all = self.manager.dict()
                self.Table = self.manager.dict()
                self.contactmap = None
                self.XL_dict = {}               

                self.model=IMP.Model()

        def compute_contact_maps_rmf3(self, rmf3, save_matrices=True):
                '''
                Compute CM for a single file.
                '''

                self.get_number_of_residues(rmf3)        
                
                self.contact_map_prob_protein_pair([rmf3], 0)
                if save_matrices:
                        self._write_matrices()

        def compute_contact_maps(self, save_matrices = True):
                '''
                Compute CM from an ensemble of models
                '''
        
                # Get rmf3 files
                self.RC = tools.ReadClustering(self.clustering_dir)
                rmfs_cluster = self.RC.get_rmfs_cluster(self.cluster)
                rmfs_cluster = random.sample(rmfs_cluster, self.number_of_models)
                
                self.rmf = rmfs_cluster[0]
                self.number_rmfs = len(rmfs_cluster)
                self.get_number_of_residues(rmfs_cluster[0])
                
                # Divide rmfs into groups
                ND = int(np.ceil(len(rmfs_cluster)/float(self.nproc)))
                rmfs_dict = {}
                for k in range(self.nproc-1):
                        rmfs_dict[k] = rmfs_cluster[(k*ND):(k*ND+ND)]
                rmfs_dict[self.nproc-1] = rmfs_cluster[((self.nproc-1)*ND):(len(rmfs_cluster))]

                # Define an output queue
                output = mp.Queue()
                
                # Setup a list of processes that we want to run
                processes = [mp.Process(target=self.contact_map_prob_protein_pair,
                                        args=(rmfs_dict[x], 0)) for x in range(self.nproc)]

                # Run processes
                for p in processes:
                        p.start()

                # Exit the completed processes
                for p in processes:
                        p.join()
                
                self._normalize_matrices()
                if save_matrices:
                        self._write_matrices()

        def get_all_indexes(self, p1):
                '''
                Get all residue indexes in bead
                '''    

                bd = IMP.atom.Fragment(p1).get_residue_indexes()
                # Check if not a fragment
                if len(bd) == 0: bd = [IMP.atom.Residue(p1).get_index()]
                return bd
        
        def get_number_of_residues(self, rmf):
                '''
                Get number of residues in each chain
                '''
                self.n_residues = {}
                model = IMP.Model()
                hier = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf)[0]
                prots = self._get_particles_at_lowest_resolution(hier)
                for p, components in prots.iteritems():
                        ri = self.get_all_indexes(components[0])[0]
                        rf = self.get_all_indexes(components[-1])[-1]
                        self.n_residues[p] = rf - ri + 1
                
	def _get_coords_array(self, particles_list):
                '''
                Get all beads coordinates and radii
                '''
		coords = []
		radii = []
		for p in particles_list:
			residue_indexes = IMP.pmi.tools.get_residue_indexes(p)
			if len(residue_indexes) !=0 :
				for res in range(min(residue_indexes), max(residue_indexes) + 1):
					d = IMP.core.XYZR(p)
					crd = np.array([d.get_x(), d.get_y(), d.get_z()])
					coords.append(crd)
					radii.append(d.get_radius())

		return np.array(coords), np.array(radii)

	def _get_contactmap_pair(self, particles_1, particles_2):
		'''
		Given two proteins, computes the contact map 
		'''
		
		coords1, radii1 = self._get_coords_array(particles_1)
		coords2, radii2  = self._get_coords_array(particles_2)
		distances = scipy.spatial.distance.cdist(coords1, coords2)
		distances = (distances - radii2).T - radii1
		#contact_map = np.where((distances>0) & (distances <= self.cutoff), 1.0, 0)
                contact_map = np.where((distances <= self.cutoff), 1.0, 0)
		return contact_map
		
	def contact_map_prob_protein_pair(self, rmfs, t):
        
		for k, rmf in enumerate(rmfs):
                        model = IMP.Model()
                        hier = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf)[0]
			prot_dictionary = self._get_particles_at_lowest_resolution(hier)
			for name1, name2 in itertools.combinations_with_replacement(prot_dictionary.keys(),2):
                                
				if name1 == name2:
                                        particles_1 = prot_dictionary[name1]
                                        particles_2 = prot_dictionary[name2]
                                        cm = self._get_contactmap_pair(particles_1,particles_2)
					if str(name1+'-'+name2) not in self.cm_all.keys():
						self.cm_all[name1+'-'+name2] = cm
					else:
						self.cm_all[name1+'-'+name2] = self.cm_all[name1+'-'+name2] + cm
				else:
                                        particles_1 = prot_dictionary[name1]
                                        particles_2 = prot_dictionary[name2]
					cm = self._get_contactmap_pair(particles_1, particles_2)
					if str(name1+'-'+name2) not in self.cm_all.keys():
						self.cm_all[name1+'-'+name2] = cm
					else:
						self.cm_all[name1+'-'+name2] = self.cm_all[name1+'-'+name2] + cm
                               
                        # Get XLs distances
                        if len(self.XL_dict)> 0:
			    self.get_XLs_distances(hier)

        def _normalize_matrices(self):
            print('Normalizing', self.number_rmfs)
            for key in self.cm_all.keys():
                    self.cm_all[key] = self.cm_all[key]/float(self.number_rmfs)
                
        def _write_matrices(self):
            '''
            Write all contact maps to files
            '''

            for key in self.cm_all.keys():
                np.savetxt(self.out_dir+'/ContMap_%s.dat'%(key),np.array(self.cm_all[key]),fmt='%s')
    

	def _get_resi_dict(self):
		self.resi_dict = {}
                model = IMP.Model()
                hier = IMP.pmi.analysis.get_hiers_from_rmf(model,0,self.rmf)[0]
                prot_dictionary = self._get_particles_at_lowest_resolution(hier)

		for prot in prot_dictionary.keys():
			resi = []
			for p in prot_dictionary[prot]:
				pp = IMP.atom.Selection(p).get_selected_particles()[0]
				r = IMP.atom.Residue(pp)
                                if 'bead' not in r.get_name():
				    resi.append([r.get_index(),str(r.get_residue_type())])
                                else:
                                    rn = r.get_name().split('-')
                                    rr = np.arange(int(rn[0]),int(rn[1].split('_')[0])+1,1)
                                    for ri in rr:
                                        resi.append([ri,'BEA'])

			self.resi_dict[prot] = resi

	def get_close_contacts(self, threshold = 0.0):
		self._get_resi_dict()

		cont_dict = {}
		for name1, name2 in itertools.combinations(self.resi_dict,2):
			contacts = []
			# Check if proteins have any contact
                        if name1+'-'+name2 in self.cm_all.keys():
                            p1 = name1
                            p2 = name2
                            mat = self.cm_all[name1+'-'+name2]  
                        elif name2+'-'+name1 in self.cm_all.keys():
                            p1 = name2
                            p2 = name1
                            mat = self.cm_all[name2+'-'+name1]
                        else:
                            raise TypeError("No contact matrix for "+ name1, name2)
 
			if np.sum(mat) > 0:
				loc = np.where(mat> threshold)
				frq = mat[loc]
				
				k = 0
				for i in np.transpose(loc):
					contacts.append(self.resi_dict[p1][i[1]]+self.resi_dict[p2][i[0]]+[int(frq[k])])
					k += 1
				# Sort array
				contacts = np.array(contacts)
				sort_contacts = contacts[contacts[:,4].astype(int).argsort()]
				np.savetxt(self.out_dir+'/contacts_%s_%s.dat'%(p1,p2),sort_contacts[::-1],fmt='%s')
				cont_dict[p1+'-'+p2] = contacts
		
	def plot_contact_maps(self,
	                      filename = 'contact_map.pdf'):
		import matplotlib as mpl
		mpl.use('Agg')
		import matplotlib.pylab as pl
		import matplotlib.gridspec as gridspec

		# Determine the number of residues per protein
		tot = 0
		nres = {}
		for p, residues in self.n_residues.iteritems():
			nres[p] = residues
			tot = tot + residues

		# Determine the proportions
		for p in nres.keys():
			nres[p] = 10*nres[p]/float(tot)

		# Create the grid
		fig = pl.figure(figsize=(8,8))
		
		nprot = len(nres.keys())
		gs = gridspec.GridSpec(nprot, nprot,
	                               width_ratios = nres.values(),
				       height_ratios = nres.values())

	
                leng = 0
		for i, p1 in enumerate(nres.keys()):
			for j, p2 in enumerate(nres.keys()):
                                print(p1, p2, nres[p1], nres[p2])
				# Get CM matrix
				if p1+'-'+p2 in self.cm_all.keys():
                                        if np.shape(self.cm_all[p1+'-'+p2])[0] == self.n_residues[p1]:
					        M = np.transpose(self.cm_all[p1+'-'+p2])
	                                else:
                                                M = np.transpose(self.cm_all[p1+'-'+p2]) 
				else:
                                        if np.shape(self.cm_all[p2+'-'+p1])[0] == self.n_residues[p2]:
					        M = self.cm_all[p2+'-'+p1]
                                        else:
                                                M = self.cm_all[p2+'-'+p1]
				ax = pl.subplot(gs[i,j])
				ax.matshow(M,cmap=pl.get_cmap('Blues'))
                                ax.set_xlim([1,np.shape(M)[1]])
                                ax.set_ylim([1,np.shape(M)[0]])
				# Plot cross-links by looping over the XL table
                                if len(self.XL_dict)> 0:
				        sXLs = { key:value for key, value in self.Table.items() if 
                                                 (((p1==key[0]) and (p2==key[2])) or ((p2==key[0]) and (p1==key[2]))) }
                                       
				        for xl in sXLs.keys():
					        pp1 = xl[0]
					        pp2 = xl[2]
					        r1 = int(xl[1])
					        r2 = int(xl[3])
                                                if len(xl)==4:
                                                        c = int(xl[4])
                                                        area = 1.5*c
                                                else:
                                                        area = 20.
					        if p1 == p2:
						        if np.min(sXLs[xl]) < self.XLs_cutoff:
							        ax.scatter(r1,r2,s=area, color='greenyellow',alpha=0.7, edgecolors='none')
							        ax.scatter(r2,r1,s=area, color='greenyellow',alpha=0.7, edgecolors='none')
						        else:
							        ax.scatter(r1,r2,s=area, color='orangered',alpha=0.7, edgecolors='none')
							        ax.scatter(r2,r1,s=area, color='orangered',alpha=0.7, edgecolors='none')
					        else:
						        if p1 == pp1 and p2 == pp2:
							        if np.min(sXLs[xl]) < self.XLs_cutoff:
								        ax.scatter(r2,r1,s=area, color='greenyellow',alpha=0.7, edgecolors='none')
							        else:
								        ax.scatter(r2,r1,s=area, color='orangered',alpha=0.7, edgecolors='none')
						        elif p2 == pp1 and p1 == pp2:
							        if np.min(sXLs[xl]) < self.XLs_cutoff:
								        ax.scatter(r1,r2,s=area, color='greenyellow',alpha=0.7, edgecolors='none')
							        else:
								        ax.scatter(r1,r2,s=area, color='orangered',alpha=0.7, edgecolors='none')

				
  				

		make_ticklabels_invisible(fig)
		# Add horizontal labels (top plots)
		for i in range(nprot):
			fig.get_axes()[i].set_title('%s'%(nres.keys()[i]),
									   fontsize=8)
		# Add vertical labels
		fig.get_axes()[0].set_ylabel('%s'%(nres.keys()[0]),
			                         fontsize=8)
		k = 1 
		for i in range(nprot,nprot*nprot,nprot):
			fig.get_axes()[i].set_ylabel('%s'%(nres.keys()[k]),
			                         fontsize=8)
			k += 1
		fig.tight_layout()
		fig.savefig(self.out_dir+'/'+filename) 


        def add_XLs(self, data_file):
            '''
            Read XLs from csv used for modeling
            '''
	    for i, line in enumerate(open(data_file)):
	        vals = line.split(',')
		if i > 0:
                        if len(vals)>4:
		                self.XL_dict[i] = [vals[0],vals[2],vals[1],vals[3],vals[5]]
                        else:
                                self.XL_dict[i] = [vals[0],vals[2],vals[1],vals[3]]

        def get_XLs_distances(self, hier):
		
		all_dist = []

		for key, value in self.XL_dict.items():
                        
			ps1 = IMP.atom.Selection(hier,
						 molecule = self.XL_dict[key][0],
						 residue_index = int(self.XL_dict[key][1]),
			                         resolution=1).get_selected_particles()
			ps2 = IMP.atom.Selection(hier,
			                         molecule=self.XL_dict[key][2],
			                         residue_index = int(self.XL_dict[key][3]),
			                         resolution=1).get_selected_particles()
			
		
			all_dist = []
			if len(ps1) == 0:
				ps1 = IMP.atom.Selection(hier,
							 molecules = self.XL_dict[key][0],
							 residue_index = int(self.XL_dict[key][1]),
							 resolution=1).get_selected_particles()
                                
			if len(ps2) == 0:
				ps2 = IMP.atom.Selection(hier,
							 molecules = self.XL_dict[key][2],
							 residue_index = int(self.XL_dict[key][3]),
							 resolution=1).get_selected_particles()
			for p1,p2 in itertools.product(ps1,ps2):
				if p1==p2:
					continue
				else:
					dist = IMP.core.get_distance(IMP.core.XYZ(p1),IMP.core.XYZ(p2))
					all_dist.append(dist)
			if len(all_dist)> 0:
			    if tuple(value) in self.Table.keys():
				self.Table[tuple(value)] = np.min([self.Table[tuple(value)],np.min(all_dist)])
			    else:
				self.Table[tuple(value)] = np.min(all_dist)
                 
        def read_contact_maps(files):
                self.cm_all = {}
                for f in files:
                        prots = f.split('ContMap_')[-1].split('.dat')[0]
                        cm = np.loadtxt(f)
                        self.cm_all[prot]

                # Get min distance for XLs
		xl_min = {}
		for i in self.Table.keys():
			xl_min[i] = np.min(self.Table[i])

        def _get_particles_at_lowest_resolution(self, hier):
                '''
                Read rmf3 file and return only coordinates of beads at
                lowest resolution
                '''
                
                particles_dict = {}
                for mol in IMP.atom.get_by_type(hier,
                                                IMP.atom.MOLECULE_TYPE):
        
                        if (len(self.selection) == 0) or (mol.get_name() in self.selection): 
                                sel = IMP.atom.Selection(mol,resolution=1)
                                particles_dict[mol.get_name()] = sel.get_selected_particles()
                
        
                return particles_dict
                         
def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)


