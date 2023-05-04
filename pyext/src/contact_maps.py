#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   Tools for contact maps analysis
"""
from __future__ import print_function
import IMP
import RMF
import IMP.algebra
import IMP.em
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.output
import IMP.rmf
import IMP.pmi.analysis

import itertools
import numpy as np
import pandas as pd
import scipy.spatial.distance
import os
import collections

import multiprocessing as mp


class CMTable(object):
    """Compute contacts maps for all models in ensemble"""

    def __init__(
        self,
        analys_dir,
        rmf_A,
        rmf_B,
        clustering_dir,
        cluster,
        selection=[],
        cutoff=16.0,
        XLs_cutoff=35.0,
        nproc=10,
    ):

        self.analys_dir = analys_dir
        self.rmf_A = rmf_A
        self.rmf_B = rmf_B
        self.clustering_dir = clustering_dir
        self.cluster = cluster
        self.selection = selection
        self.cutoff = cutoff
        self.XLs_cutoff = XLs_cutoff
        self.nproc = nproc
        self.include_XLs = False

        # Contact maps directory
        self.out_dir = os.path.join(clustering_dir, "CMs")
        if os.path.exists(self.out_dir):
            print(f"Overwriting {self.out_dir}")
        else:
            os.makedirs(self.out_dir)

        # Dictionaries for results
        self.manager = mp.Manager()
        self.cm_all = self.manager.dict()
        self.Table = self.manager.dict()
        self.contactmap = None

        self.model = IMP.Model()

    def compute_contact_maps_rmf3(self, rmf, save_matrices=True):
        """
        Compute CM from a single file.
        """

        self.get_number_of_residues(rmf)

        self.contact_map_prob_protein_pair([rmf], 0)
        if save_matrices:
            self._write_matrices()

    def compute_contact_maps(self, save_matrices=True):
        """
        Compute CM from an ensemble of models
        """

        # Number of frames in A
        model = IMP.Model()
        rh = RMF.open_rmf_file_read_only(self.rmf_A)
        n_frames_A = rh.get_number_of_frames()
        del model

        # Get frames info from clustering files
        frames_A = self.get_frames_cluster(half="A")
        frames_B = self.get_frames_cluster(half="B")
        frames_B = [f - n_frames_A for f in frames_B]

        number_of_frames = len(frames_A) + len(frames_B)

        # Divide frames into groups based on
        # number of processors
        nproc2 = int(self.nproc / 2.0)
        ND_A = int(np.ceil(len(frames_A) / float(nproc2)))
        ND_B = int(np.ceil(len(frames_B) / float(nproc2)))
        frames_A_dict = {}
        frames_B_dict = {}
        for k in range(nproc2 - 1):
            frames_A_dict[k] = frames_A[(k * ND_A):(k * ND_A + ND_A)]
            frames_A_dict[nproc2 - 1] = frames_A[
                ((nproc2 - 1) * ND_A):(len(frames_A))
            ]

        for k in range(nproc2 - 1):
            frames_B_dict[k] = frames_B[(k * ND_B):(k * ND_B + ND_B)]
            frames_B_dict[nproc2 - 1] = frames_B[
                ((nproc2 - 1) * ND_B):(len(frames_B))
            ]

        # Setup a list of processes that we want to run
        processes = [
            mp.Process(
                target=self.update_contact_maps,
                args=(self.rmf_A, frames_A_dict[x])
            )
            for x in range(nproc2)
        ] + [
            mp.Process(
                target=self.update_contact_maps,
                args=(self.rmf_B, frames_B_dict[x])
            )
            for x in range(nproc2)
        ]

        # self.contact_map_prob_protein_pair

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

        self.normalize_matrices(number_of_frames)
        if save_matrices:
            self.write_matrices()

    def get_frames_cluster(self, half="A"):
        """
        Get rmf frames for an imp-sampcon cluster
        """
        frames = []
        file_name = \
            f"{self.clustering_dir}/cluster.{self.cluster}.sample_{half}.txt"
        for line in open(
            file_name, "r"
        ):
            frames.append(int(line.strip()))
        return frames

    def get_all_indexes(self, p1):
        """
        Get all residue indexes in bead
        """

        bd = IMP.atom.Fragment(p1).get_residue_indexes()
        # Check if not a fragment
        if len(bd) == 0:
            bd = [IMP.atom.Residue(p1).get_index()]
        return bd

    def get_number_of_residues(self, rmf):
        """
        Get number of residues in each chain
        """
        self.n_residues = {}
        model = IMP.Model()
        hier = IMP.pmi.analysis.get_hiers_from_rmf(model, 0, rmf)[0]
        prots = self._get_particles_at_lowest_resolution(hier)
        for p, components in prots.items():
            ri = self.get_all_indexes(components[0])[0]
            rf = self.get_all_indexes(components[-1])[-1]
            self.n_residues[p] = rf - ri + 1

    def _get_coords_array(self, particles_list):
        """
        Get all beads coordinates and radii
        """
        coords = []
        radii = []
        for p in particles_list:
            residue_indexes = IMP.pmi.tools.get_residue_indexes(p)
            if len(residue_indexes) != 0:
                for res in range(min(residue_indexes),
                                 max(residue_indexes) + 1):
                    d = IMP.core.XYZR(p)
                    crd = np.array([d.get_x(), d.get_y(), d.get_z()])
                    coords.append(crd)
                    radii.append(d.get_radius())

        return np.array(coords), np.array(radii)

    def _get_contactmap_pair(self, particles_1, particles_2):
        """
        Given two proteins, computes the contact map
        """

        coords1, radii1 = self._get_coords_array(particles_1)
        coords2, radii2 = self._get_coords_array(particles_2)
        distances = scipy.spatial.distance.cdist(coords1, coords2)
        distances = (distances - radii2).T - radii1
        contact_map = np.where((distances <= self.cutoff), 1.0, 0)
        return contact_map

    def update_contact_maps(self, rmf_file, frames):
        """
        Iterate over protein pairs and compute distance
        matrices
        """
        model = IMP.Model()
        rh = RMF.open_rmf_file_read_only(rmf_file)
        hier = IMP.rmf.create_hierarchies(rh, model)[0]

        # Iterate over frames
        for fr in frames:
            IMP.rmf.load_frame(rh, RMF.FrameID(fr))
            prot_dictionary = self._get_particles_at_lowest_resolution(hier)
            for name1, name2 in itertools.combinations_with_replacement(
                prot_dictionary.keys(), 2
            ):
                if name1 == name2:
                    particles_1 = prot_dictionary[name1]
                    particles_2 = prot_dictionary[name2]
                    cm = self._get_contactmap_pair(particles_1, particles_2)
                    if f"{name1}-{name2}" not in self.cm_all.keys():
                        self.cm_all[f"{name1}-{name2}"] = cm
                    else:
                        self.cm_all[f"{name1}-{name2}"] += cm
                else:
                    particles_1 = prot_dictionary[name1]
                    particles_2 = prot_dictionary[name2]
                    cm = self._get_contactmap_pair(particles_1, particles_2)
                    if f"{name1}-{name2}" not in self.cm_all.keys():
                        self.cm_all[f"{name1}-{name2}"] = cm
                    else:
                        self.cm_all[f"{name1}-{name2}"] += cm
            # Get XLs distances
            if self.include_XLs:
                self.get_XLs_distances(hier)

    def normalize_matrices(self, number_of_frames):
        """
        Normalize matrices by
        number of models
        """
        print("Normalizing ...")
        for key in self.cm_all.keys():
            self.cm_all[key] = self.cm_all[key] / float(number_of_frames)

    def write_matrices(self):
        """
        Write all contact maps to files
        """

        for key in self.cm_all.keys():
            np.savetxt(
                os.path.join(self.out_dir, f"ContMap_{key}.dat"),
                np.array(self.cm_all[key]),
                fmt="%s",
            )

    def _get_resi_dict(self):
        self.resi_dict = {}
        model = IMP.Model()
        hier = IMP.pmi.analysis.get_hiers_from_rmf(model, 0, self.rmf)[0]
        prot_dictionary = self._get_particles_at_lowest_resolution(hier)

        for prot in prot_dictionary.keys():
            resi = []
            for p in prot_dictionary[prot]:
                pp = IMP.atom.Selection(p).get_selected_particles()[0]
                r = IMP.atom.Residue(pp)
                if "bead" not in r.get_name():
                    resi.append([r.get_index(), str(r.get_residue_type())])
                else:
                    rn = r.get_name().split("-")
                    rr = np.arange(int(rn[0]), int(rn[1].split("_")[0]) + 1, 1)
                    for ri in rr:
                        resi.append([ri, "BEA"])

            self.resi_dict[prot] = resi

    def get_close_contacts(self, threshold=0.2):
        self._get_resi_dict()

        cont_dict = {}
        for name1, name2 in itertools.combinations(self.resi_dict, 2):
            contacts = []
            # Check if proteins have any contact
            if f"{name1}-{name2}" in self.cm_all.keys():
                p1 = name1
                p2 = name2
                mat = self.cm_all[f"{name1}-{name2}"]
            elif f"{name1}-{name2}" in self.cm_all.keys():
                p1 = name2
                p2 = name1
                mat = self.cm_all[f"{name1}-{name2}"]
            else:
                raise TypeError("No contact matrix for " + name1, name2)

            if np.sum(mat) > 0:
                loc = np.where(mat > threshold)
                frq = mat[loc]

                k = 0
                for i in np.transpose(loc):
                    contacts.append(
                        self.resi_dict[p1][i[1]] +
                        self.resi_dict[p2][i[0]] +
                        [frq[k]]
                    )
                    k += 1
                # Sort array
                contacts = np.array(contacts)
                sort_contacts = \
                    contacts[contacts[:, 4].astype(float).argsort()]
                np.savetxt(
                    os.path.join(self.out_dir, f"contacts_{p1}_{p2}.dat"),
                    sort_contacts[::-1],
                    fmt="%s",
                )
                cont_dict[f"{p1}-{p2}"] = contacts

    def plot_contact_maps(self, scaling_factor=10.,
                          filename="contact_map.pdf"):
        """
        Plot all protein pairs contact maps
        scaling_factor is used to set the area of the XLs dots
        """

        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pylab as pl
        import matplotlib.gridspec as gridspec

        self.get_number_of_residues(self.rmf_A)

        # Determine the number of residues per protein
        tot = 0
        nres = {}
        for p, residues in self.n_residues.items():
            nres[p] = residues
            tot = tot + residues

        # Determine the proportions
        for p in nres.keys():
            nres[p] = 10 * nres[p] / float(tot)

        # Create the grid
        fig = pl.figure(figsize=(8, 8))

        nprot = len(nres.keys())
        gs = gridspec.GridSpec(
            nprot,
            nprot,
            width_ratios=nres.values(),
            height_ratios=nres.values()
        )

        for i, p1 in enumerate(nres.keys()):
            for j, p2 in enumerate(nres.keys()):
                # Get CM matrix
                if f"{p1}-{p2}" in self.cm_all.keys():
                    M = np.transpose(self.cm_all[f"{p1}-{p2}"])
                else:
                    M = self.cm_all[f"{p2}-{p1}"]
                ax = pl.subplot(gs[i, j])
                ax.matshow(M, cmap=pl.get_cmap("Blues"))
                ax.set_xlim([1, np.shape(M)[1]])
                ax.set_ylim([1, np.shape(M)[0]])

                # Plot cross-links by looping over the XL table
                if self.include_XLs:
                    sXLs = {
                        key: value
                        for key, value in self.Table.items()
                        if (p1.split(".")[0] == key[0] and
                            p2.split(".")[0] == key[1])
                        or (p2.split(".")[0] == key[0] and
                            p1.split(".")[0] == key[1])
                    }

                    for xl in sXLs.keys():
                        pp1 = xl[0]
                        pp2 = xl[1]
                        r1 = int(xl[2])
                        r2 = int(xl[3])
                        if len(xl) > 4:
                            c = float(xl[4])
                            area = scaling_factor * c
                        else:
                            area = scaling_factor
                        if p1 == p2:
                            if np.min(sXLs[xl]) < self.XLs_cutoff:
                                ax.scatter(
                                    r1,
                                    r2,
                                    s=area,
                                    color="limegreen",
                                    alpha=0.7,
                                    edgecolors="none",
                                )
                                ax.scatter(
                                    r2,
                                    r1,
                                    s=area,
                                    color="limegreen",
                                    alpha=0.7,
                                    edgecolors="none",
                                )
                            else:
                                ax.scatter(
                                    r1,
                                    r2,
                                    s=area,
                                    color="orangered",
                                    alpha=0.7,
                                    edgecolors="none",
                                )
                                ax.scatter(
                                    r2,
                                    r1,
                                    s=area,
                                    color="orangered",
                                    alpha=0.7,
                                    edgecolors="none",
                                )
                        else:
                            if (p1.split(".")[0] == pp1
                                    and p2.split(".")[0] == pp2):
                                if np.min(sXLs[xl]) < self.XLs_cutoff:
                                    ax.scatter(
                                        r2,
                                        r1,
                                        s=area,
                                        color="limegreen",
                                        alpha=0.7,
                                        edgecolors="none",
                                    )
                                else:
                                    ax.scatter(
                                        r2,
                                        r1,
                                        s=area,
                                        color="orangered",
                                        alpha=0.7,
                                        edgecolors="none",
                                    )
                            elif (p2.split(".")[0] == pp1
                                  and p1.split(".")[0] == pp2):
                                if np.min(sXLs[xl]) < self.XLs_cutoff:
                                    ax.scatter(
                                        r1,
                                        r2,
                                        s=area,
                                        color="limegreen",
                                        alpha=0.7,
                                        edgecolors="none",
                                    )
                                else:
                                    ax.scatter(
                                        r1,
                                        r2,
                                        s=area,
                                        color="orangered",
                                        alpha=0.7,
                                        edgecolors="none",
                                    )

        make_ticklabels_invisible(fig)

        # Add horizontal labels (top plots)
        for i in range(nprot):
            fig.get_axes()[i].set_title(
                "%s" % (list(nres.keys())[i]), fontsize=8
            )
        # Add vertical labels
        fig.get_axes()[0].set_ylabel(
            "%s" % (list(nres.keys())[0]), fontsize=8
        )
        k = 1
        for i in range(nprot, nprot * nprot, nprot):
            fig.get_axes()[i].set_ylabel(
                "%s" % (list(nres.keys())[k]), fontsize=8
            )
            k += 1
        fig.tight_layout()
        fig.savefig(os.path.join(self.out_dir, filename))

    def plot_contact_maps_subunits(self):
        """
        Function to plot distances between beads
        /residues between all combinations of
        subunits
        """
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pylab as pl
        import matplotlib.pyplot as plt

        plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.labelbottom"] = True
        plt.rcParams["xtick.top"] = plt.rcParams["xtick.labeltop"] = True
        # Determine the number of residues per protein
        tot = 0
        nres = collections.OrderedDict()
        for p, residues in self.n_residues.items():
            nres[p] = residues
            tot = tot + residues

        # Determine the proportions
        for p in nres.keys():
            nres[p] = 10 * nres[p] / float(tot)

        # Create the grid
        fig = pl.figure(figsize=(8, 8))
        combinations = list(itertools.combinations(nres.keys(), 2))
        for i, p in enumerate(combinations):
            p1 = p[0]
            p2 = p[1]
            fig = pl.figure(figsize=(8, 8))
            filename = f"Contact_map_{p1}-{p2}.pdf"
            if f"{p1}-{p2}" in self.cm_all.keys():
                if np.shape(self.cm_all[f"{p1}-{p2}"])[0] == \
                   self.n_residues[p1]:
                    M = np.transpose(self.cm_all[f"{p1}-{p2}"])
                else:
                    M = np.transpose(self.cm_all[f"{p1}-{p2}"])
            else:
                if np.shape(self.cm_all[f"{p2}-{p1}"])[0] == \
                   self.n_residues[p2]:
                    M = self.cm_all[f"{p2}-{p1}"]
                else:
                    M = self.cm_all[f"{p2}-{p1}"]
            ax = fig.add_subplot(111)
            cax = ax.matshow(M, cmap=pl.get_cmap("RdYlBu"), vmin=0, vmax=100)
            ax.set_xlim([1, np.shape(M)[1]])
            ax.set_ylim([1, np.shape(M)[0]])
            ax.tick_params(
                axis="x",
                bottom=True,
                top=False,
                labelbottom=True,
                labeltop=False,
                labelsize=20,
            )
            ax.tick_params(axis="y", labelsize=20)
            ax.set_xlabel(p2, fontsize=20)
            ax.set_ylabel(p1, fontsize=20)
            fig.colorbar(cax, ax=ax, ticks=[0, 20, 40, 60, 80, 100])
            # cmap.ax.tick_params(labelsize=14)
            fig.tight_layout()
            fig.savefig(
                os.path.join(self.out_dir, filename), bbox_inches="tight"
            )
            plt.close("all")
            plt.close(fig)
            plt.clf()

    def add_XLs_data(
        self, data_file, keys=["Protein1", "Protein2", "Residue1", "Residue2"]
    ):
        """
        Read XLs from csv used for modeling
        """
        D = pd.read_csv(data_file)
        D_XLs = D[keys]

        # Order columns
        prot_columns = [
            c for c in D_XLs.columns.values
            if isinstance(D_XLs.iloc[0][c], str)
        ]
        other = [c for c in D_XLs.columns.values if c not in prot_columns]

        D_XLs = D[prot_columns + other]

        # Rename columns
        D_XLs = D_XLs.rename(
            columns={D_XLs.columns[0]: "Protein A",
                     D_XLs.columns[1]: "Protein B",
                     D_XLs.columns[2]: "Residue A",
                     D_XLs.columns[3]: "Residue B"})
        self.D_XLs = D_XLs

        self.include_XLs = True

    def get_XLs_distances(self, hier):

        all_dist = []

        for i, row in self.D_XLs.iterrows():

            ps1 = IMP.atom.Selection(
                hier,
                molecule=row["Protein A"],
                residue_index=int(row["Residue A"]),
                resolution=1,
            ).get_selected_particles()
            ps2 = IMP.atom.Selection(
                hier,
                molecule=row["Protein B"],
                residue_index=int(row["Residue B"]),
                resolution=1,
            ).get_selected_particles()

            if len(ps1) > 0 and len(ps2) > 0:
                for p1, p2 in itertools.product(ps1, ps2):
                    if p1 == p2:
                        continue
                    else:
                        dist = IMP.core.get_distance(IMP.core.XYZ(p1),
                                                     IMP.core.XYZ(p2))
                        all_dist.append(dist)
                if len(all_dist) > 0:
                    value = row[:].values
                    if tuple(value) in self.Table.keys():
                        self.Table[tuple(value)] = np.min(
                            [self.Table[tuple(value)], np.min(all_dist)]
                        )
                    else:
                        self.Table[tuple(value)] = np.min(all_dist)

    def read_contact_maps(self, files):
        """
        Read contact map from text file
        """
        self.cm_all = {}
        for f in files:
            prot = f.split("ContMap_")[-1].split(".dat")[0]
            cm = np.loadtxt(f)
            self.cm_all[prot] = cm

        # Get min distance for XLs
        xl_min = {}
        for i in self.Table.keys():
            xl_min[i] = np.min(self.Table[i])

    def _get_particles_at_lowest_resolution(self, hier):
        """
        Read rmf3 file and return only coordinates of beads at
        lowest resolution
        """

        particles_dict = {}
        for mol in IMP.atom.get_by_type(hier, IMP.atom.MOLECULE_TYPE):

            copy = IMP.atom.Copy(mol).get_copy_index()
            if len(self.selection) == 0 or mol.get_name() in self.selection:
                sel = IMP.atom.Selection(mol, resolution=1)
                particles_dict[
                    f"{mol.get_name()}.{copy}"
                ] = sel.get_selected_particles()

        return particles_dict


def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)
