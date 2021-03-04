#!/usr/bin/env python

'''
Compute different distance metrics to
compare structures in a ensemble
'''

import os
import random
import math
import itertools
import numpy as np
import multiprocessing as mp
import pandas as pd

import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import RMF

import tools


class get_distance_metrics(object):
    def __init__(self, clustering_dir, cluster, number_of_models=200,
                 selection=[], align_to=[], write_all_values=False):

        self.clustering_dir = clustering_dir
        self.cluster = cluster
        self.selection = selection
        self.align_to = align_to
        self.number_of_models = number_of_models
        self.write_all_values = write_all_values

        self.manager = mp.Manager()
        self.all_coords = self.manager.dict()
        self.lock = self.manager.Lock()
        self.all_versus_all = 0
        self.all_versus_centroid = 0

        self.RC = tools.ReadClustering(self.clustering_dir)

        rmfs_cluster = self.RC.get_rmfs_cluster(self.cluster)

        if len(rmfs_cluster) > self.number_of_models:
            rmfs_cluster = random.sample(rmfs_cluster, self.number_of_models)
        self.structures = rmfs_cluster

        self.get_rbs_from_centroid()

        self._read_all_coords_RBs(rmfs_cluster)

        for k, v in self.rb_components.items():
            if len(v) > 1:
                name = ''
                for vv in v:
                    name += str(vv[0])+' '+str(vv[1])+'/'
                name = name[0:-1]
            else:
                name = str(v[0][0])+' '+str(v[0][1])
            print(k, name)

    def compute_DRMSD_all_versus_all(self):
        i = 0
        ntot = len(list(itertools.combinations(self.all_coords.values(), 2)))
        for particles1, particles2 in itertools.combinations(
                self.all_coords.values(), 2):
            self.compute_DRMSD_pairs_mp(particles1, particles2)
            i += 1
            if i % 100 == 0:
                print(str(i) + ' of ' + str(ntot))

        self.all_versus_all = 1
        self.write_DRMSD_output()

    def compute_DRMSD_all_versus_centroid(self):
        rmf3 = os.path.join(self.clustering_dir, 'cluster.'+str(self.cluster),
                            'cluster_center_model.rmf3')

        particles1 = self._get_RBs_particle_coordinates(rmf3)

        ntot = len(self.all_coords.values())

        for i, particles2 in enumerate(self.all_coords.values()):
            self.compute_DRMSD_pairs_mp(particles1, particles2)
            if i % 100 == 0:
                print(str(i) + ' of ' + str(ntot))

        self.all_versus_centroid = 1
        self.write_DRMSD_output()

    def write_DRMSD_output(self):

        if self.all_versus_all == 1:
            file_matrix = 'DRMSD_mean_pairs_all_versus_all.csv'
            file_index = 'RBs_indexes_DRMSD_calculation_all_versus_all.dat'
            if self.write_all_values:
                file_values = 'DRMSD_pairs_all_versus_all_'

        elif self.all_versus_centroid == 1:
            file_matrix = 'DRMSD_mean_pairs_all_versus_centroid.csv'
            file_index \
                = 'RBs_indexes_DRMSD_calculation_all_versus_centroid.dat'
            if self.write_all_values:
                file_values = 'DRMSD_pairs_all_versus_centroid_'

        else:
            raise ValueError('ERROR: Need to compute DRMSD')

        # Create matrix for plotting
        M = pd.DataFrame(0, columns=self.rb_components.keys(),
                         index=self.rb_components.keys())
        for pair, values in self.DRMSD_pairs.items():
            M.loc[pair[0], pair[1]] = np.mean(values)
            M.loc[pair[1], pair[0]] = M.loc[pair[0], pair[1]]
            if self.write_all_values:
                np.savetxt(file_values+str(pair[0])+'_'+str(pair[1])+'.csv',
                           np.array(values))
        print(M)

        M.to_csv(file_matrix)

        out = open(file_index, 'w')
        for k, v in self.rb_components.items():
            out.write(str(k)+'\t'+str(v).strip('[]')+'\n')
        out.close()

    def compute_RMSD_all_versus_centroid(self):

        self.RMSD_pairs = self.manager.dict()
        for i, j in itertools.combinations_with_replacement(
                self.rb_components.keys(), 2):
            self.RMSD_pairs[(i, j)] = []

        rmf3 = self.clustering_dir + '/cluster.' + \
            str(self.cluster) + '/cluster_center_model.rmf3'

        particles1, reference1 = self._get_RBs_particle_coordinates(rmf3)
        ntot = len(self.structures)
        print(ntot)

        rmfs_cluster = list(zip(range(len(self.structures)), self.structures))
        ND = int(np.ceil(len(rmfs_cluster)/float(self.nproc)))

        rmfs_dict = {}
        for k in range(self.nproc - 1):
            rmfs_dict[k] = rmfs_cluster[(k*ND):(k*ND+ND)]
        rmfs_dict[self.nproc-1] \
            = rmfs_cluster[((self.nproc-1)*ND):(len(rmfs_cluster))]

        if self.nproc > 1:
            processes = [mp.Process(target=self._compute_RMSD_mp,
                                    args=(particles1, reference1,
                                          rmfs_dict[x]))
                         for x in range(self.nproc)]
            for p in processes:
                p.start()

            for p in processes:
                p.join()

        elif self.nproc == 1:
            self._compute_RMSD_mp(particles1, reference1, rmfs_dict[0])

        self.all_versus_centroid = 1
        self.write_RMSD_output()

    def _compute_RMSD_mp(self, particles1, reference1, rmfs_dict_sub):

        for i, structure in rmfs_dict_sub:
            particles2, reference2 \
                = self._get_RBs_particle_coordinates(structure)
            self.compute_RMSD(particles1, reference1, particles2, reference2)

    def compute_RMSD(self, particles1, reference1, particles2, reference2):

        import pyRMSD.RMSDCalculator
        calculator_name = "QCP_SERIAL_CALCULATOR"

        if set(particles1.keys()) != set(particles2.keys()):
            raise ValueError('ERROR: rmfs have different RBs')

        else:

            for rb1, rb2 in itertools.combinations_with_replacement(
                    particles1.keys(), 2):
                cc = np.array([particles1[rb1] + particles1[rb2],
                               particles2[rb1] + particles2[rb2]])

                cr = np.array([reference1['ref'],
                               reference2['ref']])

                calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
                    calculator_name, fittingCoordsets=cr,
                    calculationCoordsets=cc)

                rmsd = calculator.pairwise(
                    0, 1, get_superposed_coordinates=False)
                with self.lock:
                    self.RMSD_pairs[(rb1, rb2)] += [rmsd]

    def write_RMSD_output(self):
        out_dir = os.path.join(self.clustering_dir, 'RMSDs')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        file_matrix = os.path.join(
            out_dir, 'RMSD_mean_pairs_all_versus_centroid.csv')
        file_index = os.path.join(
            out_dir, 'RBs_indexes_RMSD_calculation_all_versus_centroid.dat')
        file_values = 'RMSD_pairs_all_versus_centroid_'

        M = pd.DataFrame(0, columns=self.rb_components.keys(),
                         index=self.rb_components.keys())
        for pair, values in self.RMSD_pairs.items():
            M.loc[pair[0], pair[1]] = np.mean(values)
            M.loc[pair[1], pair[0]] = M.loc[pair[0], pair[1]]
            file_out = os.path.join(
                out_dir, file_values+str(pair[0])+'_'+str(pair[1])+'.csv')

            np.savetxt(file_out, np.array(values))
        print(M)

        M.to_csv(file_matrix)

        out = open(file_index, 'w')
        for k, v in self.rb_components.items():
            out.write(str(k)+'\t'+str(v).strip('[]')+'\n')
        out.close()

    def _read_all_coords_RBs(self, rmfs_cluster):
        '''
        Multiprocessing reading of all
        RMFs
        '''

        self.nproc = 10
        rmfs_cluster = zip(range(len(rmfs_cluster)), rmfs_cluster)

        # Divide rmfs into groups
        ND = int(np.ceil(len(rmfs_cluster)/float(self.nproc)))
        rmfs_dict = {}
        for k in range(self.nproc-1):
            rmfs_dict[k] = rmfs_cluster[(k*ND):(k*ND+ND)]
        rmfs_dict[self.nproc-1] = rmfs_cluster[
            ((self.nproc-1)*ND):(len(rmfs_cluster))]

        # Setup a list of processes that we want to run
        processes = [mp.Process(target=self._read_coords_RBs,
                                args=(rmfs_dict[x], 0))
                     for x in range(self.nproc)]

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

    def _read_coords_RBs(self, numbered_rmfs, a):
        for i, file in numbered_rmfs:
            particles = self._get_RBs_particle_coordinates(file)
            self.all_coords[i] = particles

    def get_rbs_from_centroid(self):
        rmf3 = os.path.join(self.clustering_dir,
                            'cluster.'+str(self.cluster),
                            'cluster_center_model.rmf3')
        self._read_representation(rmf3)

    def _read_representation(self, rmf3):
        m = IMP.Model()
        rh_ref = RMF.open_rmf_file_read_only(rmf3)
        h_ref = IMP.rmf.create_hierarchies(rh_ref, m)[0]
        IMP.rmf.load_frame(rh_ref, RMF.FrameID(0))

        self.components = []
        self.struct_components = {}
        self.flex_components = {}
        rb_components = {}

        n_struct, n_flex = 0, 0
        for state in h_ref.get_children():
            for component in state.get_children():
                if len(self.selection) == 0 \
                        or component.get_name() in self.selection:
                    self.components.append(component.get_name())
                    for fragment in component.get_children():
                        leaves = IMP.atom.get_leaves(fragment)
                        residues = IMP.atom.Fragment(
                            fragment).get_residue_indexes()
                        residue_range = str(residues[0])+'-'+str(residues[-1])
                        n_res = residues[-1] - residues[0]
                        if not IMP.core.RigidMember.get_is_setup(leaves[0]):
                            n_flex += n_res
                            if component.get_name() in \
                                    self.flex_components.keys():
                                self.flex_components[component.get_name()] \
                                    += ', '+residue_range
                            else:
                                self.flex_components[component.get_name()] \
                                    = residue_range
                        else:
                            n_struct += n_res
                            for res in fragment.get_children():
                                p = res.get_children()[0]
                                rb = IMP.core.RigidMember(p).get_rigid_body()
                                if rb in rb_components.keys():
                                    rb_components[rb].append(
                                        (component.get_name(), residue_range))
                                else:
                                    rb_components[rb] = [(component.get_name(),
                                                          residue_range)]

        self.struc_coverage = 100.0 * n_struct / (n_struct + n_flex)

        self.rb_components = {
            i: v for i, (k, v) in enumerate(rb_components.items())}

        print('Coverage', self.struc_coverage)
        print('Number of rigid-bodies', len(self.rb_components))

        del m

    def _get_RBs_particle_coordinates(self, rmf3):
        model = IMP.Model()
        h = IMP.pmi.analysis.get_hiers_from_rmf(model, 0, rmf3)[0]

        particles_dict = {}
        for rb, sel in self.rb_components.items():
            particles = []
            for comp in self.rb_components[rb]:
                if len(self.selection) == 0 or (comp[0] in self.selection):
                    r1 = int(comp[1].split('-')[0])
                    r2 = int(comp[1].split('-')[1])
                    sel = IMP.atom.Selection(h,
                                             molecule=comp[0],
                                             residue_indexes=range(r1, r2),
                                             resolution=1)
                    particles += sel.get_selected_particles()

            coords = [np.array(IMP.core.XYZ(p).get_coordinates())
                      for p in particles]
            particles_dict[rb] = coords

        if len(self.align_to) > 0:
            reference_dict = {}
            for ref in self.align_to:
                sel = IMP.atom.Selection(
                    h, molecule=ref[0],
                    residue_indexes=range(int(ref[1]), int(ref[2])),
                    resolution=1)
                particles += sel.get_selected_particles()
            coords = [np.array(IMP.core.XYZ(p).get_coordinates())
                      for p in particles]
            reference_dict['ref'] = coords

            del model
            del h

            return particles_dict, reference_dict

        else:
            del model
            del h

            return particles_dict

    def compute_DRMSD_pairs(self, particles1, particles2):

        if set(particles1.keys()) != set(particles2.keys()):
            raise ValueError('ERROR: rmfs have different RBs')

        else:
            for rb1, rb2 in itertools.combinations(particles1.keys(), 2):
                p1_rb1 = particles1[rb1]
                p1_rb2 = particles1[rb2]
                p2_rb1 = particles2[rb1]
                p2_rb2 = particles2[rb2]
                if len(p1_rb1) == len(p2_rb1) and len(p1_rb2) == len(p2_rb2):
                    drmsd = self.get_DRMSD(p1_rb1, p1_rb2, p2_rb1, p2_rb2)
                    self.DRMSD_pairs[(rb1, rb2)] += [drmsd]
                else:
                    ValueError('ERROR: RBs have different number of particles')

    def compute_DRMSD_pairs_mp(self, particles1, particles2):

        self.DRMSD_pairs = self.manager.dict()
        for i, j in itertools.combinations(self.rb_components.keys(), 2):
            self.DRMSD_pairs[(i, j)] = []

        pairs = []
        if set(particles1.keys()) != set(particles2.keys()):
            raise ValueError('ERROR: rmfs have different RBs')

        else:
            for rb1, rb2 in itertools.combinations(particles1.keys(), 2):
                p1_rb1 = particles1[rb1]
                p1_rb2 = particles1[rb2]
                p2_rb1 = particles2[rb1]
                p2_rb2 = particles2[rb2]
                pairs.append((rb1, rb2, p1_rb1, p1_rb2, p2_rb1, p2_rb2))

        # Setup a list of processes that we want to run
        processes = [mp.Process(target=self.get_DRMSD,
                                args=(p)) for p in pairs]

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

    def get_DRMSD(self, rb1, rb2, p1_rb1, p1_rb2, p2_rb1, p2_rb2, only_RBs=1):

        import scipy.spatial.distance

        dist1 = scipy.spatial.distance.cdist(p1_rb1, p1_rb2)
        dist2 = scipy.spatial.distance.cdist(p2_rb1, p2_rb2)

        dist1 = dist1.flatten()
        dist2 = dist2.flatten()

        drmsd = math.sqrt(np.sum((dist1-dist2)**2)/len(dist1))

        self.DRMSD_pairs[(rb1, rb2)] += [drmsd]

    def _get_structured_beads(self, rmf3):
        '''
        Read rmf3 file and return only coordinates of beads at
        one residue-per-bead
        '''

        model = IMP.Model()
        h = IMP.pmi.analysis.get_hiers_from_rmf(model, 0, rmf3)[0]

        particles_dict = {}
        for mol in IMP.atom.get_by_type(h, IMP.atom.MOLECULE_TYPE):
            sel = IMP.atom.Selection(mol, resolution=1)
            particles = [p for p in sel.get_selected_particles()
                         if not IMP.atom.Fragment.get_is_setup(p)]
            coords, radii = [], []
            for p in particles:
                coords.append(np.array(IMP.core.XYZR(p).get_coordinates()))
                radii.append(IMP.core.XYZR(p).get_radius())

            particles_dict[mol.get_name()] = (np.array(coords), radii)
        return particles_dict

    def adjacent_values(self, vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value

    def set_axis_style(self, ax, labels):
        ax.get_xaxis().set_tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(np.arange(1, len(labels) + 1))
        ax.set_xticklabels(labels, rotation=45, fontsize=10)
        ax.set_xlim(0.25, len(labels) + 0.75)
        ax.set_xlabel('Sample name')
