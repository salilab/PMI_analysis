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
from operator import itemgetter
from copy import deepcopy
from math import log,sqr
import itertools
import numpy as np
import scipy.spatial.distance
import os


class CMTable(object):
	""" Compute contacts maps for all models in ensemble """
	def __init__(self,
                    cutoff = 20.0,
                    out_dir ='./'):

		self.contactmap = None
		self.cutoff = cutoff
                self.out_dir = out_dir
                self.Table = {}
                self.XL_dict = {}               
 
                if os.path.exists(self.out_dir):
                    print('Overwriting '+ self.out_dir)
                else:
                    os.makedirs(self.out_dir)
    
        def get_all_indexes(self, p1):
            '''
            Get all residue indexes in bead
            '''    

            bd = IMP.atom.Fragment(p1).get_residue_indexes()
            # Check if not a fragment
            if len(bd) == 0: bd = [IMP.atom.Residue(p1).get_index()]
            return bd

        def get_number_of_residues(self, model, rmf):
            '''
            Get number of residues in each chain
            '''
            self.n_residues = {}
            h = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf)[0]
            prots = get_particles_at_resolution_one(h)
            for p, components in prots.iteritems():
                ri = self.get_all_indexes(components[0])[0]
                rf = self.get_all_indexes(components[-1])[-1]
                self.n_residues[p] = rf - ri + 1
        
	def _get_coords_array(self, particles_list):
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
					#if name not in self.index_dictionary:
					#	self.index_dictionary[name] = [resindex]
					#else:
					#self.index_dictionary[name].append(resindex)
                    #resindex += 1

		return np.array(coords), np.array(radii)

	def _get_contactmap_pair(self, name1,name2):
		'''
		Given two proteins, computes the contact map 
		'''
		
		particles_1 = self.prot_dictionary[name1]
		particles_2 = self.prot_dictionary[name2]
		coords1, radii1 = self._get_coords_array(particles_1)
		coords2, radii2  = self._get_coords_array(particles_2)
		distances = scipy.spatial.distance.cdist(coords1, coords2)
		distances = (distances - radii2).T - radii1
		contact_map = np.where((distances>0) & (distances <= self.cutoff), 1.0, 0)

		return contact_map
		
	def contact_map_prob_protein_pair(self, model, rmfs, save_matrices = True):
		self.cm_all = {}
		k = 0
		for k, rmf in enumerate(rmfs):
			prots = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf)[0]
			self.prot_dictionary = get_particles_at_resolution_one(prots)
	
			for name1, name2 in itertools.combinations_with_replacement(self.prot_dictionary.keys(),2):
				if name1 == name2:
					cm = self._get_contactmap_pair(name1,name2)
					if str(name1+'-'+name2) not in self.cm_all.keys():
						self.cm_all[name1+'-'+name2] = cm
					else:
						self.cm_all[name1+'-'+name2] = self.cm_all[name1+'-'+name2] + cm
				else:
					cm = self._get_contactmap_pair(name1,name2)
					if str(name1+'-'+name2) not in self.cm_all.keys():
						self.cm_all[name1+'-'+name2] = cm
					else:
						self.cm_all[name1+'-'+name2] = self.cm_all[name1+'-'+name2] + cm
                        # Get XLs distances
                        if len(self.XL_dict)> 0:
			    self.get_dists(prots)
                
		# Write all contact maps to files
		if save_matrices == True:
			for key in self.cm_all.keys():
				np.savetxt(os.path.join(
                                        self.out_dir,'ContMap_%s.dat'%(key)),np.array(self.cm_all[key])/float(k),fmt='%s')

	def _get_resi_dict(self):
		self.resi_dict = {}
		for prot in self.prot_dictionary.keys():
			resi = []
			for p in self.prot_dictionary[prot]:
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

	def get_close_contacts(self, threshold = 0.5):
		self._get_resi_dict()

		cont_dict = {}
		for name1, name2 in itertools.combinations_with_replacement(self.resi_dict,2):
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
				#print('Number of contacts between:', name1, name2, 'is', int(np.sum(self.cm_all[name1+'-'+name2])))
				k = 0
				for i in np.transpose(loc):
					contacts.append(self.resi_dict[p1][i[1]]+self.resi_dict[p2][i[0]]+[int(frq[k])])
					k += 1
				# Sort array
				contacts = np.array(contacts)
				sort_contacts = contacts[contacts[:,4].astype(int).argsort()]
				np.savetxt(os.path.join(
                                        self.out_dir,'contacts_%s_%s.dat'%(p1,p2)),sort_contacts[::-1],fmt='%s')
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

		i=0
                leng = 0
		for p1 in nres.keys():
			for p2 in nres.keys():
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
                                print(p1, p2, np.shape(M))
				ax = pl.subplot(gs[i])
                                #ax.set_aspect('equal')
				ax.matshow(M,cmap=pl.get_cmap('Blues'))
                                ax.set_xlim([1,np.shape(M)[1]])
                                ax.set_ylim([1,np.shape(M)[0]])
				# Plot cross-links by looping over the XL table
				sXLs = { key:value for key, value in self.Table.items() if 
                                        (((p1==key[0]) and (p2==key[2])) or ((p2==key[0]) and (p1==key[2]))) }
                                leng += len(sXLs)
                                print('leng', len(sXLs), len(self.Table))
				for xl in sXLs.keys():
					pp1 = xl[0]
					pp2 = xl[2]
					r1 = int(xl[1])
					r2 = int(xl[3])
					if p1 == p2:
						if np.min(sXLs[xl]) < 35.0:
							ax.scatter(r1,r2, color='greenyellow',alpha=0.7)
							ax.scatter(r2,r1, color='greenyellow',alpha=0.7)
						else:
							ax.scatter(r1,r2, color='orangered',alpha=0.7)
							ax.scatter(r2,r1, color='orangered',alpha=0.7)
					else:
						if p1 == pp1 and p2 == pp2:
							if np.min(sXLs[xl]) < 35.0:
								ax.scatter(r2,r1, color='greenyellow',alpha=0.7)
							else:
								ax.scatter(r2,r1, color='orangered',alpha=0.7)
						elif p2 == pp1 and p1 == pp2:
							if np.min(sXLs[xl]) < 35.0:
								ax.scatter(r1,r2, color='greenyellow',alpha=0.7)
							else:
								ax.scatter(r1,r2, color='orangered',alpha=0.7)

				
  				i += 1

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
		fig.savefig(os.path.join(self.out_dir,filename)) 


	def coordinates_to_map(self,prots):
		coords = []
		radii = []
		namelist = []

		particles_dictionary = get_particles_at_resolution_one(prots)

		resindex = 0
		self.index_dictionary={}

		self._check_ambiguity(prots)
		
		for name in particles_dictionary.keys():
			print('prot name: ', name)
			# To do: test for ambiguity and expand dictionary
			residue_index = []
			for p in particles_dictionary[name]:
				residue_indexes = IMP.pmi.tools.get_residue_indexes(p)
				if len(residue_indexes) !=0 :
					for res in range(min(residue_indexes), max(residue_indexes) + 1):
						d = IMP.core.XYZR(p)
						crd = np.array([d.get_x(), d.get_y(), d.get_z()])
						coords.append(crd)
						radii.append(d.get_radius())
						if name not in self.index_dictionary:
							self.index_dictionary[name] = [resindex]
						else:
							self.index_dictionary[name].append(resindex)
                        resindex += 1
	
		coords = np.array(coords)
		radii = np.array(radii)
		import scipy.spatial.distance
		distances = scipy.spatial.distance.cdist(coords, coords)
		distances = (distances - radii).T - radii
		distances = np.where((distances>0) & (distances <=16), 1.0, 0)
		if self.contactmap is None:
			self.contactmap = np.zeros((len(coords), len(coords)))
		self.contactmap += distances
		print('CM: ', prots, sum(self.contactmap))
	
		return np.array(self.contactmap)

        def get_XLdict(self,data_file):
            '''
            Read XLs from csv used for modeling
            '''
	    for i, line in enumerate(open(data_file)):
	        vals = line.split(',')
		if i > 0:
		    self.XL_dict[i] = [vals[0],vals[2],vals[1],vals[3]]


        def get_dists(self, prot):
		
		all_dist = []
	
		for keys in self.XL_dict.keys():
                        
			self.prot_dictionary = get_particles_at_resolution_one(prot)
			ps1 = IMP.atom.Selection(prot,
						 molecule = self.XL_dict[keys][0],
						 residue_index = int(self.XL_dict[keys][1]),
			                         resolution=1).get_selected_particles()
			ps2 = IMP.atom.Selection(prot,
			                         molecule=self.XL_dict[keys][2],
			                         residue_index = int(self.XL_dict[keys][3]),
			                         #atom_type=IMP.atom.AT_CA,
			                         resolution=1).get_selected_particles()
			
		
			all_dist = []
			if len(ps1) == 0:
				ps1 = IMP.atom.Selection(prot,
							 molecules = self.XL_dict[keys][0],
							 residue_index = int(self.XL_dict[keys][1]),
							 resolution=1).get_selected_particles()
                                
			if len(ps2) == 0:
				ps2 = IMP.atom.Selection(prot,
							 molecules = self.XL_dict[keys][2],
							 residue_index = int(self.XL_dict[keys][3]),
							 resolution=1).get_selected_particles()
			for p1,p2 in itertools.product(ps1,ps2):
				if p1==p2:
					continue
				else:
					dist = IMP.core.get_distance(IMP.core.XYZ(p1),IMP.core.XYZ(p2))
					all_dist.append(dist)
                        #print(keys, min(all_dist))
			if len(all_dist)> 0:
				if self.XL_dict[keys][0]+'-'+self.XL_dict[keys][1]+'-'+self.XL_dict[keys][2]+'-'+self.XL_dict[keys][3] in self.Table.keys():
					self.Table[(self.XL_dict[keys][0],+self.XL_dict[keys][1],self.XL_dict[keys][2],self.XL_dict[keys][3])].append(min(all_dist))
				else:
					self.Table[(self.XL_dict[keys][0],self.XL_dict[keys][1],self.XL_dict[keys][2],self.XL_dict[keys][3])] = [min(all_dist)]									
        
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

def get_particles_at_resolution_one(prot):
    """
    Get particles at res 1, or any beads, based on the name.
    No Representation is needed. This is mainly used when the hierarchy
    is read from an RMF file.
    @return a dictionary of component names and their particles
    \note If the root node is named "System" or is a "State", do proper selection.
    """
    particle_dict = {}

    # attempt to give good results for PMI2
    if IMP.pmi.get_is_canonical(prot):
        for mol in IMP.atom.get_by_type(prot,IMP.atom.MOLECULE_TYPE):
            sel = IMP.atom.Selection(mol,resolution=1)
            particle_dict[mol.get_name()] = sel.get_selected_particles()
    else:
        allparticles = []
        for c in prot.get_children():
            name = c.get_name()
            particle_dict[name] = IMP.atom.get_leaves(c)
            for s in c.get_children():
                if "_Res:1" in s.get_name() and "_Res:10" not in s.get_name():
                    allparticles += IMP.atom.get_leaves(s)
                if "Beads" in s.get_name():
                    allparticles += IMP.atom.get_leaves(s)

        particle_align = []
        for name in particle_dict:
            particle_dict[name] = IMP.pmi.tools.sort_by_residues(
                list(set(particle_dict[name]) & set(allparticles)))
    return particle_dict

def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)

