#!/usr/bin/env python

"""
Tools to analyze PMI MC runs
"""

from __future__ import division

import os
import math
import glob
import shutil
import random
import itertools
import subprocess
import pandas as pd
import numpy as np
import multiprocessing as mp
from equilibration import detectEquilibration
import hdbscan

import IMP
import IMP.rmf
import RMF

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pylab as pl  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402

mpl.rcParams.update({"font.size": 10})


color_palette = [
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
    "#ccb974",
    "#64b5cd",
    "#4c72b0",
    "#55a868",
    "#c44e52",
    "#8172b2",
]


class AnalysisTrajectories(object):
    def __init__(
        self,
        out_dirs,
        dir_name="run_",
        analysis_dir="analysis",
        detect_equilibration=True,
        burn_in_fraction=0.02,
        nskip=20,
        nproc=1,
        number_models_out=29999,
        plot_fmt="pdf",
    ):
        """
        Analyze the ensemble of models obtained after structural sampling by
        filtering out the bad scoring models, and clustering the good scoring
        models according to one or multiple score terms, thereby obtaining
        a set of representative good scoring models that can be used to
        subsequently compute sampling precision and model precision [using the
        sampcon module].

        @param out_dirs: list of directories containing the trajectory (rmf
        files and stat files) for each independent run.

        @param dir_name: prefix for all directories containing structural
        sampling runs. Usually these are named as <prefix>1, <prefix>2, etc.
        Default is "run_"

        @analysis_dir: name of output directory where all analysis results will
        be stored. Default is "analysis".

        @detect_equilibration: if True, the equilibration period is rigorously
        determined from the trajectory of all scores using statistical
        inefficiency measures. Default is True.

        @burn_in_fraction: a <burn_in_fraction> fraction of frames are always
        discarded from the beginning of the trajectory, when calculating
        statistics. Default is 0.02 (i.e. 2 % of the total trajectory size).

        @nskip: skips every <nskip> frames when calculating statistics.
        Default is 200.

        @nproc: number of processors for parallel execution (uses the
        python multiprocess module). Default is 1.

        @plot_fmt: file extensions for all plots created. Default is pdf.
        """

        self.out_dirs = [os.path.abspath(d) for d in out_dirs]
        self.dir_name = dir_name
        self.analysis_dir = os.path.abspath(analysis_dir)
        self.detect_equilibration = detect_equilibration
        self.nproc = nproc
        self.nskip = nskip
        self.burn_in_frac = burn_in_fraction
        self.number_models_out = number_models_out

        # check if plot_fmt is supported
        supported_fmts = list(pl.gcf().canvas.get_supported_filetypes().keys())
        if plot_fmt not in supported_fmts:
            raise KeyError(
                "plot_fmt not found in supported matplotlib " "file types:",
                supported_fmts,
            )
        self.plot_fmt = plot_fmt

        self.restraint_names = {}
        self.all_score_fields = []
        self.rerun = False

        # report if equilibration detection has been requested
        if not self.detect_equilibration:
            print("Not running equilibration check.")

        # report burn-in
        if self.burn_in_frac > 0:
            print(
                "Considering %1.2f fraction of frames from the beginning "
                "of each independent run as the burn-in phase and "
                "discarding them" % self.burn_in_frac
            )

        # Create analysis dir if missing
        if not os.path.isdir(self.analysis_dir):
            os.mkdir(self.analysis_dir)

        # For multiprocessing
        self.manager = mp.Manager()
        self.S_all = self.manager.dict()
        self.S_info_all = self.manager.dict()
        self.S_dist_all = self.manager.dict()
        self.XLs_nuis = self.manager.dict()

        # Global sampling parameters
        self.Sampling = self.manager.dict()
        s_vals = ["Number_of_replicas", "N_equilibrated", "N_total"]

        for v in s_vals:
            self.Sampling[v] = []
        self.Validation = self.manager.dict()

        # Define with restraints to analyze
        self.Connectivity_restraint = False
        self.Excluded_volume_restraint = False
        self.XLs_restraint = False
        self.XLs_restraint_nuisances = True
        self.Multiple_XLs_restraints = False
        self.atomic_XLs_restraint = False
        self.atomic_XLs_restraint_nuisances = True
        self.Multiple_atomic_XLs_restraints = False

        self.EM_restraint = False
        self.Distance_restraint = False
        self.Binding_restraint = False
        self.Occams_restraint = False
        self.Occams_restraint_nuisances = False
        self.Occams_positional_restraint = False
        self.Occams_positional_nuisances = False
        self.pEMAP_restraint = False
        self.pEMAP_restraint_new = False
        self.DOPE_restraint = False
        self.MembraneExclusion_restraint = False
        self.MembraneSurfaceLocation_restraint = False
        self.score_only_restraint = False

        # Other handles
        self.restraint_handles = [
            "MonteCarlo_Nframe", "rmf_frame_index", "Total_Score"]
        self.restraint_names = {"MonteCarlo_Nframe": "MC_frame"}
        self.restraint_names["rmf_frame_index"] = "rmf_frame_index"
        self.restraint_names["Total_Score"] = "Total_Score"
        self.distance_handles = []
        self.info_handles = []

        # By default, add restraints of same type
        self.sum_Connectivity_restraint = True
        self.sum_Excluded_volume_restraint = True
        self.sum_Binding_restraint = True
        self.sum_Distance_restraint = True
        self.sum_XLs_restraint = True
        self.sum_Occams_positional_restraint = True
        self.sum_atomic_XLs_restraint = True
        self.sum_DOPE_restraint = True
        self.sum_MembraneExclusion_restraint = True
        self.sum_MembraneSurfaceLocation_restraint = True
        self.sum_EM_restraint = False
        self.sum_score_only_restraint = {}

        self.Multiple_psi_values = False

        # Separate trajectories into two halves
        self.dir_halfA = np.sort(self.out_dirs)[::2]
        self.dir_halfB = np.sort(self.out_dirs)[1::2]

    def set_analyze_XLs_restraint(
        self,
        get_nuisances=True,
        Multiple_XLs_restraints=False,
        ambiguous_XLs_restraint=False,
        XLs_cutoffs={"DSSO": 30.0},
    ):

        self.restraint_handles.append(
            "CrossLinkingMassSpectrometryRestraint_Data_Score"
        )
        self.restraint_handles.append(
            "CrossLinkingMassSpectrometryRestraint_PriorPsi_Score"
        )
        self.restraint_names[
            "CrossLinkingMassSpectrometryRestraint_Data_Score"] = "XLs"
        self.restraint_names[
            "CrossLinkingMassSpectrometryRestraint_PriorPsi_Score"
        ] = "XLs_psi"

        self.distance_handles.append(
            "CrossLinkingMassSpectrometryRestraint_Distance_")
        self.info_handles.append("CrossLinkingMassSpectrometryRestraint_Psi_")

        # self.XLs_restraint = True
        self.Multiple_XLs_restraints = Multiple_XLs_restraints
        if self.Multiple_XLs_restraints:
            self.sum_XLs_restraint = False
        self.ambiguous_XLs_restraint = ambiguous_XLs_restraint

        self.XLs_cutoffs = XLs_cutoffs
        self.XLs_restraint = True

    def set_analyze_atomic_XLs_restraint(
        self,
        get_nuisances=True,
        Multiple_atomic_XLs_restraints=False,
        atomic_XLs_cutoffs={"DSSO": 30.0},
    ):
        self.atomic_XLs_restraint = True
        self.select_atomic_XLs_satisfaction = True
        if get_nuisances:
            self.atomic_XLs_restraint_nuisances = True
        if Multiple_atomic_XLs_restraints:
            self.Multiple_atomic_XLs_restraints = True
            self.sum_atomic_XLs_restraint = False
        self.atomic_XLs_cutoffs = atomic_XLs_cutoffs
        self.ambiguous_XLs_restraint = False

    def set_analyze_Connectivity_restraint(self):
        self.restraint_handles.append("ConnectivityRestraint")
        self.restraint_names["ConnectivityRestraint"] = "CR"

    def set_analyze_Excluded_volume_restraint(self):
        self.restraint_handles.append("ExcludedVolumeSphere")
        self.restraint_names["ExcludedVolumeSphere"] = "EV"

    def set_analyze_EM_restraint(self):
        self.restraint_handles.append("GaussianEMRestraint")
        self.restraint_names["GaussianEMRestraint"] = "EM3D"
        self.EM_restraint = True

    def set_analyze_Distance_restraint(self):
        self.restraint_handles.append("DistanceRestraint_Score")
        self.restraint_names["DistanceRestraint_Score"] = "DR"
        self.Distance_restraint = True

    def set_analyze_Binding_restraint(self):
        self.restraint_handles.append("ResidueBindingRestraint_score")
        self.restraint_names["ResidueBindingRestraint_score"] = "BR"
        self.Binding_restraint = True

    def set_analyze_Occams_restraint(self):
        self.restraint_handles.append("OccamsRestraint_Score")
        self.restraint_handles.append("OccamsRestraint_psi_Score")
        self.restraint_handles.append("OccamsRestraint_sigma_Score")
        self.restraint_names["OccamsRestraint_Score"] = "Struct_Equiv"
        self.restraint_names["OccamsRestraint_psi_Score"] = "Struct_Equiv_psi"
        self.restraint_names["OccamsRestraint_sigma_Score"] \
            = "Struct_Equiv_sigma"

        # Other relevant info
        self.info_handles.append("OccamsRestraint_satisfied")
        self.info_handles.append("OccamsRestraint_sigma")
        self.info_handles.append("OccamsRestraint_psi")
        self.Occams_restraint = True

    def set_analyze_Occams_positional_restraint(self):
        self.restraint_handles.append("OccamsPositionalRestraint_Score")
        self.restraint_names["OccamsPositionalRestraint_Score"] = "OccPos"

        self.Occams_positional_restraint = True
        # self.Occams_positional_nuisances = True

    def set_analyze_pEMAP_restraint(self):
        self.restraint_handles.append("SimplifiedPEMAP_data_Score")
        self.restraint_names["SimplifiedPEMAP_data_Score"] = "pEMap"
        # self.pEMAP_restraint = True

    def set_analyze_pEMAP_restraint_new(self):
        self.restraint_handles.append("pEMapRestraint_Score")
        self.restraint_names["pEMapRestraint_Score"] = "pEMap"
        self.info_handles.append("pEMapRestraint_satisfaction")
        self.info_handles.append("pEMapRestraint_sigma")
        self.pEMAP_restraint_new = True

    def set_analyze_DOPE_restraint(self):
        self.restraint_handles.append("DOPE_Restraint_score")
        self.restraint_names["DOPE_Restraint_score"] = "DOPE"
        # self.DOPE_restraint = True

    def set_analyze_MembraneExclusion_restraint(self):
        self.restraint_handles.append("MembraneExclusionRestraint")
        self.restraint_names["MembraneExclusionRestraint"] = "MEX"
        # self.MembraneExclusion_restraint = True

    def set_analyze_MembraneSurfaceLocation_restraint(self):
        self.restraint_handles.append("MembraneSurfaceLocation")
        self.restraint_names["MembraneSurfaceLocation"] = "MSL"
        # self.MembraneSurfaceLocation_restraint = True

    def set_analyze_score_only_restraint(self, handle, short_name,
                                         do_sum=True):
        self.restraint_handles.append(handle)
        self.restraint_names[handle] = short_name
        self.score_only_restraint = True
        self.sum_score_only_restraint[short_name] = do_sum

    def set_select_by_Total_score(self, score_cutoff):
        self.select_Total_score = True
        self.cutoff_Total_score = score_cutoff

    def set_select_by_EM_score(self, score_cutoff):
        self.select_EM_score = True
        self.cutoff_EM_score = score_cutoff

    def get_score_fields(self, stat2_dict):
        score_fields = []
        score_names = []

        for restraint in self.restraint_handles:
            RES = {
                stat2_dict[k]: k
                for k in stat2_dict.keys()
                if (restraint in stat2_dict[k])
                and ("Acceptance" not in stat2_dict[k])
                and ("StepSize" not in stat2_dict[k])
            }

            if len(RES) == 1:
                if restraint not in [
                    "MonteCarlo_Nframe",
                    "rmf_frame_index",
                    "Total_Score",
                ]:
                    score_names.append(
                        self.restraint_names[restraint] + "_sum")
                else:
                    score_names.append(self.restraint_names[restraint])
                score_fields += RES.values()
            else:
                for k, v in RES.items():
                    score_names.append(
                        self.restraint_names[restraint]
                        + "_"
                        + k.split(restraint + "_")[-1]
                    )
                    score_fields += [v]

        return score_fields, score_names

    def get_distance_fields(self, stat2_dict):
        """Get inter-beads distances fields"""
        dist_dict = {}

        for handle in self.distance_handles:
            dist = {
                stat2_dict[k]: k
                for k in stat2_dict.keys()
                if (handle in stat2_dict[k] and "Score" not in stat2_dict[k])
            }
            dist_dict.update(dist)
        dist_names = list(dist_dict.keys())
        dist_names.sort()
        dist_fields = [dist_dict[k] for k in dist_names]

        return list(dist_names), dist_fields

    def get_info_fields(self, stat2_dict):
        """Get other relevent fields
        (ex. percent satisfaction, nuisance values)"""
        info_dict = {}
        for handle in self.info_handles:
            info = {
                stat2_dict[k]: k
                for k in stat2_dict.keys()
                if (
                    handle in stat2_dict[k]
                    and "MonteCarlo_" not in stat2_dict[k]
                    and "Score" not in stat2_dict[k]
                )
            }
            info_dict.update(info)
        info_names = list(info_dict.keys())
        info_names.sort()
        info_fields = [info_dict[k] for k in info_names]

        return list(info_names), info_fields

    def read_DB(self, db_file):
        """Read database"""
        DB = {}
        i = 0
        with open(db_file, "r") as of:
            db_lines = of.readlines()
        for line in db_lines:
            vals = line.split()
            DB[str(i)] = vals
            i += 1
        return DB

    def get_keys(self, stat_file):
        """Get all keys in stat file"""

        with open(stat_file, "r") as of:
            stat_file_lines = of.readlines()
        for line in stat_file_lines:
            d = eval(line)
            klist = list(d.keys())
            # check if it is a stat2 file
            if "STAT2HEADER" in klist:
                import operator

                for k in klist:
                    if "STAT2HEADER" in str(k):
                        del d[k]
                stat2_dict = d
                # get the list of keys sorted by value
                kkeys = [
                    k[0] for k in sorted(stat2_dict.items(),
                                         key=operator.itemgetter(1))
                ]
                klist = [
                    k[1] for k in sorted(stat2_dict.items(),
                                         key=operator.itemgetter(1))
                ]
                invstat2_dict = {}
                for k in kkeys:
                    invstat2_dict.update({stat2_dict[k]: k})
            else:
                klist.sort()
            break

        return stat2_dict

    def read_stat_files(self):
        """Multiprocessor reading of stat files"""

        self.Sampling["Number_of_runs"] = len(self.out_dirs)

        # Split directories
        ND = int(np.ceil(len(self.out_dirs) / float(self.nproc)))
        out_dirs_dict = {}
        for k in range(self.nproc - 1):
            out_dirs_dict[k] = list(self.out_dirs[(k * ND): (k * ND + ND)])
        out_dirs_dict[self.nproc - 1] = list(
            self.out_dirs[((self.nproc - 1) * ND): (len(self.out_dirs))]
        )

        # Setup a list of processes that we want to run
        processes = [
            mp.Process(target=self.read_traj_info, args=((out_dirs_dict[x],)))
            for x in range(self.nproc)
        ]

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

    def read_stats_detailed(self, traj, stat_files):
        """
        Detailed reading of stats files that includes
        the rmf in which the frame is.
        To be used when using rmf_slice
        """

        S_scores = []
        S_dist = []
        P_info = []

        for sf in stat_files:
            # Read header
            stat2_dict = self.get_keys(sf)
            # Get fields to extract
            query_rmf_file = self.get_field_id(stat2_dict, "rmf_file")

            # Compile fields to read
            score_fields, score_names = self.get_score_fields(stat2_dict)
            dist_names, dist_fields = self.get_distance_fields(stat2_dict)
            info_names, info_fields = self.get_info_fields(stat2_dict)

            line_number = 0
            with open(sf, "r") as of:
                sf_lines = of.readlines()

            for line in sf_lines:
                line_number += 1
                try:
                    d = eval(line)
                except:  # noqa: E722
                    print(
                        "# Warning: skipped line number "
                        + str(line_number)
                        + " not a valid line"
                    )
                    break

                if line_number > 1:
                    frmf = [d[field] for field in query_rmf_file][0]
                    s0 = [float(d[field]) for field in score_fields] \
                        + [traj, frmf]
                    S_scores.append(s0)
                    if len(dist_fields) > 0:
                        d0 = [s0[0]] + [float(d[field])
                                        for field in dist_fields]
                        S_dist.append(d0)
                    if len(info_fields) > 0:
                        p0 = [s0[0]] + [float(d[field])
                                        for field in info_fields]
                        P_info.append(p0)

        # Sort based on frame
        S_scores.sort(key=lambda x: float(x[0]))

        # Convert into pandas DF
        column_names = [x for x in score_names] + ["traj", "rmf3_file"]
        DF = pd.DataFrame(S_scores, columns=column_names)

        # If some restraints need to be added
        if self.XLs_restraint and self.sum_XLs_restraint:
            XLs_sum = pd.Series(self.add_restraint_type(DF, "XLs_"))
            DF = DF.assign(XLs_sum=XLs_sum.values)

        if self.atomic_XLs_restraint and self.sum_atomic_XLs_restraint:
            atomic_XLs_sum = pd.Series(
                self.add_restraint_type(DF, "atomic_XLs_"))
            DF = DF.assign(atomic_XLs_sum=atomic_XLs_sum.values)

        all_EV = [res for res in score_names if "EV_" in res]
        if len(all_EV) > 1:
            EV_sum = pd.Series(self.add_restraint_type(DF, "EV_"))
            DF = DF.assign(EV_sum=EV_sum.values)

        all_CR = [res for res in score_names if "CR_" in res]
        if len(all_CR) > 1:
            CR_sum = pd.Series(self.add_restraint_type(DF, "CR_"))
            DF = DF.assign(CR_sum=CR_sum.values)

        all_BR = [res for res in score_names if "BR_" in res]
        if len(all_BR) > 1:
            BR_sum = pd.Series(self.add_restraint_type(DF, "BR_"))
            DF = DF.assign(BR_sum=BR_sum.values)

        if self.Distance_restraint and self.sum_Distance_restraint:
            DR_sum = pd.Series(self.add_restraint_type(DF, "DR_"))
            DF = DF.assign(DR_sum=DR_sum.values)
        if (self.MembraneExclusion_restraint
                and self.sum_MembraneExclusion_restraint):
            MEX_sum = pd.Series(self.add_restraint_type(DF, "MEX_"))
            DF = DF.assign(MEX_sum=MEX_sum.values)

        if (self.MembraneExclusion_restraint
                and self.sum_MembraneExclusion_restraint):
            MSL_sum = pd.Series(self.add_restraint_type(DF, "MSL_"))
            DF = DF.assign(MSL_sum=MSL_sum.values)

        # if some score only type custom restraints defined by the user
        # need to be summed. Because of the do_sum flag available while
        # calling the set_analyze_score_only_restraint() method on the
        # user defined restraint, they don't need to separately handled
        # in read_traj_info. If you don't want the sum to be analyzed, just
        # don't ever compute it by setting do_sum = False.
        if self.score_only_restraint:
            for key, do_sum in self.sum_score_only_restraint.items():
                prefix = key + "_"
                all_restraint = [res for res in score_names if prefix in res]
                if len(all_restraint) > 1 and do_sum:
                    this_restraint_sum = pd.Series(
                        self.add_restraint_type(DF, prefix))
                    sum_dict = {prefix + "sum": this_restraint_sum.values}
                    DF = DF.assign(**sum_dict)

        # Get distance fields
        if len(dist_names) > 0:
            S_dist = np.array(S_dist)
            S_dist = S_dist[S_dist[:, 0].argsort()]
            # Convert in DF
            DF_dXLs = pd.DataFrame(S_dist, columns=["MC_frame"] + dist_names)

        if len(info_names) > 0:
            P_info = np.array(P_info)
            P_info = P_info[P_info[:, 0].argsort()]
            DF_info = pd.DataFrame(P_info, columns=["MC_frame"] + info_names)

        # return DF, DF_XLs, P_satif, frames_dic
        if dist_fields and info_fields:
            return DF, DF_dXLs, DF_info
        elif dist_fields and not info_fields:
            return DF, DF_dXLs, None
        elif not dist_fields and info_fields:
            return DF, None, DF_info
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

        if isinstance(out_dirs_sel, str):
            out_dirs_sel = [out_dirs_sel]

        for out in out_dirs_sel:
            if self.dir_name in out:
                traj = [x for x in out.split("/") if self.dir_name in x][0]
                traj_number = int(traj.split(self.dir_name)[1])
            else:
                traj = 0
                traj_number = 0
            stat_files = sorted(glob.glob(os.path.join(out, "stat.*.out")))
            self.Sampling["Number_of_replicas"] += [len(stat_files)]

            # Read all stat files of trajectory
            S_tot_scores, S_dist, S_info = self.read_stats_detailed(traj,
                                                                    stat_files)

            n_frames = len(S_tot_scores)
            burn_in = int(self.burn_in_frac * n_frames)
            print(
                "The mean score, min score, and n frames are: ",
                traj_number,
                np.mean(S_tot_scores["Total_Score"].iloc[burn_in:]),
                np.min(S_tot_scores["Total_Score"].iloc[burn_in:]),
                n_frames,
            )
            self.Sampling["N_total"] = self.Sampling["N_total"] + [n_frames]

            # Selection of just sums (default)
            sel_entries = ["Total_Score"] + [
                v for v in S_tot_scores.columns.values if "sum" in v
            ]

            # If specified by user, can look at individual contributions
            if (self.Connectivity_restraint
                    and not self.sum_Connectivity_restraint):
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "CR_" in v and "sum" not in v
                ]
            if (
                self.Excluded_volume_restraint
                and not self.sum_Excluded_volume_restraint
            ):
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "EV_" in v and "sum" not in v
                ]
            if self.Binding_restraint and not self.sum_Binding_restraint:
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "BR_" in v and "sum" not in v
                ]
            if self.Distance_restraint and not self.sum_Distance_restraint:
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "DR_" in v and "sum" not in v
                ]
            if self.XLs_restraint and not self.sum_XLs_restraint:
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "XLs_" in v and "sum" not in v
                ]
            if self.atomic_XLs_restraint and not self.sum_atomic_XLs_restraint:
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "atomic_XLs_" in v and "sum" not in v
                ]
            if self.DOPE_restraint and not self.sum_DOPE_restraint:
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "DOPE_" in v and "sum" not in v
                ]
            if self.EM_restraint and not self.sum_EM_restraint:
                sel_entries += [
                    v
                    for v in S_tot_scores.columns.values
                    if "EM3D_" in v and "sigma" not in v
                ]

            # Also add nuisances parameter
            sel_entries += [
                v for v in S_tot_scores.columns.values
                if "Psi" in v and "sum" not in v
            ]

            # Detect equilibration time if requested
            ts_eq = []
            if self.detect_equilibration:
                for r in sel_entries:
                    try:
                        [t, g, N] = detectEquilibration(
                            np.array(S_tot_scores[r].loc[burn_in:]),
                            nskip=self.nskip
                        )
                        ts_eq.append(t)
                    except (ValueError, TypeError):
                        ts_eq.append(0)
                print("Trajectory, ts_eqs: ", traj, ts_eq)
            else:
                ts_eq = [0] * len(sel_entries)

            ts_max = np.max(ts_eq) + burn_in
            self.Sampling["N_equilibrated"] += [n_frames - ts_max]

            # Plot the scores and restraint satisfaction
            file_out = "plot_scores_%s.%s" % (traj_number, self.plot_fmt)
            self.plot_scores_restraints(
                S_tot_scores[["MC_frame"] + sel_entries],
                ts_eq, burn_in, file_out
            )

            if self.pEMAP_restraint_new:
                file_out_pemap = "plot_pEMAP_%s.%s" % (traj_number,
                                                       self.plot_fmt)
                self.plot_pEMAP_satisfaction(S_info, file_out_pemap)

            if self.Occams_restraint:
                file_out_occams = "plot_Occams_satisfaction_%s.%s" % (
                    traj_number,
                    self.plot_fmt,
                )
                self.plot_Occams_satisfaction(S_info, file_out_occams)

            # Check how many XLs are satisfied
            if self.XLs_restraint:
                S_info = self.analyze_trajectory_XLs(
                    S_dist,
                    S_info,
                    atomic_XLs=False,
                    traj_number=traj_number,
                    ts_max=ts_max,
                )

            if self.atomic_XLs_restraint:
                S_info = self.analyze_trajectory_XLs(
                    S_dist,
                    S_info,
                    atomic_XLs=True,
                    traj_number=traj_number,
                    ts_max=ts_max,
                )

            # Add half info
            if out in self.dir_halfA:
                S_tot_scores = S_tot_scores.assign(
                    half=pd.Series(
                        ["A"] * len(S_tot_scores), index=S_tot_scores.index
                    ).values
                )
            elif out in self.dir_halfB:
                S_tot_scores = S_tot_scores.assign(
                    half=pd.Series(
                        ["B"] * len(S_tot_scores), index=S_tot_scores.index
                    ).values
                )
            else:
                S_tot_scores = S_tot_scores.assign(
                    half=pd.Series(
                        [0] * len(S_tot_scores), index=S_tot_scores.index
                    ).values
                )

            # Collect distances and nuisances information

            self.S_all[out] = S_tot_scores[ts_max:]
            try:
                self.S_info_all[out] = S_info[ts_max:]
            except:  # noqa: E722
                print("No S_info")
            if self.XLs_restraint:
                self.S_dist_all[out] = S_dist[ts_max:]
            if self.atomic_XLs_restraint:
                self.S_dist_all[out] = S_dist[ts_max:]

    def get_field_id(self, dict, val):
        """
        For single field, get number of fields in stat file
        """
        return [k for k in dict.keys() if dict[k] == val]

    def plot_scores_restraints(self, selected_scores, ts_eq, burn_in,
                               file_out):
        """
        For each trajectory plot all restraint scores
        """
        n_bins = 20
        ts_max = np.max(ts_eq)
        n_res = len(selected_scores.columns.values) - 1

        fig, ax = pl.subplots(figsize=(2.0 * n_res, 4.0), nrows=2, ncols=n_res)
        axes = ax.flatten()
        for i, c in enumerate(selected_scores.columns.values[1:]):
            axes[i].plot(
                selected_scores["MC_frame"].loc[burn_in::10],
                selected_scores[c].loc[burn_in::10],
                color="b",
                alpha=0.5,
            )
            axes[i].axvline(ts_eq[i], color="grey")
            axes[i].set_title(c, fontsize=14)
            axes[i].set_xlabel("Step", fontsize=12)
            if i == 0:
                axes[i].set_ylabel("Score (a.u.)", fontsize=12)

        for i, c in enumerate(selected_scores.columns.values[1:]):
            axes[i + n_res].hist(
                selected_scores[c].loc[ts_eq[i]:: 10],
                n_bins,
                histtype="step",
                fill=False,
                color="orangered",
                alpha=0.9,
            )
            axes[i + n_res].hist(
                selected_scores[c].loc[ts_max::10],
                n_bins,
                histtype="step",
                fill=False,
                color="gold",
                alpha=0.9,
            )
            axes[i + n_res].set_xlabel("Score (a.u.)", fontsize=12)
            if i == 0:
                axes[i + n_res].set_ylabel("Density", fontsize=12)

        pl.tight_layout(pad=0.5, w_pad=0.1, h_pad=2.0)
        fig.savefig(os.path.join(self.analysis_dir, file_out))
        pl.close()

    def write_models_info(self):
        """
        Write info of all models after equilibration
        """

        for k, T in self.S_all.items():
            kk = k.split(self.dir_name)[-1].split("/")[0]
            T.to_csv(
                os.path.join(
                    os.path.join(self.analysis_dir,
                                 "scores_info_" + str(kk) + ".csv")
                )
            )

        for k in self.S_dist_all.keys():
            T = self.S_dist_all[k]
            kk = k.split(self.dir_name)[-1].split("/")[0]
            T.to_csv(
                os.path.join(
                    os.path.join(self.analysis_dir,
                                 "XLs_dist_info_" + str(kk) + ".csv")
                )
            )

        for k in self.S_info_all.keys():
            T = self.S_info_all[k]
            kk = k.split(self.dir_name)[-1].split("/")[0]
            T.to_csv(
                os.path.join(
                    os.path.join(self.analysis_dir,
                                 "other_info_" + str(kk) + ".csv")
                )
            )

    def read_models_info(self, XLs_cutoffs=None):
        """
        Read info of all models after equilibration
        """

        self.rerun = True
        if XLs_cutoffs:
            self.XLs_cutoffs = XLs_cutoffs

        # Score files
        info_files = glob.glob(os.path.join(self.analysis_dir,
                                            "scores_info_*.csv"))
        for f in info_files:
            k = f.split("all_info_")[-1].split(".csv")[0]
            df = pd.read_csv(f)
            self.S_all[k] = df

        # XLs files
        xls_files = glob.glob(os.path.join(self.analysis_dir,
                                           "XLs_dist_info_*.csv"))
        if len(xls_files) > 0:
            self.XLs_restraint = True
            self.ambiguous_XLs_restraint = False
            for f in xls_files:
                k = f.split("XLs_info_")[-1].split(".csv")[0]
                df = pd.read_csv(f)
                self.S_dist_all[k] = df

            # Check for ambiguity
            k0 = list(self.S_dist_all.keys())[0]
            XLs_names = self.S_dist_all[k0].columns.values
            self.ambiguous_XLs_dict = self.check_XLs_ambiguity(XLs_names)
        else:
            print("No files with XLs info found")

    def hdbscan_clustering(
        self, selected_scores, min_cluster_size=150, min_samples=5, skip=1
    ):
        """
        DO HDBSCAN clustering for selected restraint and/or nuisance parameters
        """

        all_dfs = [self.S_all[dd] for dd in np.sort(self.S_all.keys())]
        S_comb = pd.concat(all_dfs)

        # Print all available fields before checking if field exists.
        print("All available fields: ", S_comb.columns.values)
        S_comb_sel = S_comb[selected_scores].iloc[::skip]
        S_comb_all = S_comb.iloc[::skip]
        print("Fields selected for HDBSCAN clustering: ",
              S_comb_sel.columns.values)

        hdbsc = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size, min_samples=min_samples
        ).fit(S_comb_sel)

        # Add clusters labels
        S_comb_sel = S_comb_sel.assign(
            cluster=pd.Series(hdbsc.labels_, index=S_comb_sel.index).values
        )
        S_comb_all = S_comb_all.assign(
            cluster=pd.Series(hdbsc.labels_, index=S_comb_all.index).values
        )

        # Add cluster labels also to XLs info if available
        if self.XLs_restraint is True:

            all_dist_dfs = [
                self.S_dist_all[dd] for dd in np.sort(self.S_dist_all.keys())
            ]
            S_comb_dist_clustering = pd.concat(
                all_dist_dfs, sort=False).iloc[::skip]

            S_comb_dist_clustering = S_comb_dist_clustering.assign(
                cluster=pd.Series(
                    hdbsc.labels_, index=S_comb_dist_clustering.index
                ).values
            )
            S_comb_dist_clustering.to_csv(
                os.path.join(self.analysis_dir, "XLs_clustering_info.csv"),
                index=False)
            self.S_comb_dist_clustering = S_comb_dist_clustering

        print(
            "Number of unique clusters: ",
            len(np.unique(hdbsc.labels_)),
            np.unique(hdbsc.labels_),
        )

        # Write and plot info from clustering
        self.plot_hdbscan_clustering(S_comb_sel, selected_scores)
        self.write_hdbscan_clustering(S_comb_all)
        self.plot_hdbscan_runs_info(S_comb_all)

        S_comb_sel = S_comb_sel.assign(
            half=pd.Series(S_comb_all["half"], index=S_comb_sel.index).values
        )
        self.write_summary_hdbscan_clustering(
            S_comb_all[["Total_Score"] + selected_scores + ["half", "cluster"]]
        )

    def write_summary_hdbscan_clustering(self, S_comb_all):
        """
        Write clustering summary information
        (i.e. cluster number, average scores, number of models)
        """

        aggregation = {
            v: "mean" for v in S_comb_all.columns.values
            if v not in ["half", "cluster"]
        }
        aggregation.update(
            {"cluster": lambda x: len(x), "half": lambda x: len(x[x == "A"])}
        )

        S_clusters = S_comb_all.groupby("cluster").agg(aggregation)

        S_clusters.rename(columns={"cluster": "N_models", "half": "N_A"},
                          inplace=True)
        S_clusters["N_B"] = S_clusters["N_models"] - S_clusters["N_A"]
        S_clusters = S_clusters.sort_values("Total_Score")

        # Selection criteria: lowest total score with at least 20-80 split
        for i, row in S_clusters.iterrows():
            if 1.0 - (row["N_models"] - row["N_A"]) / row["N_models"] >= 0.2:
                self.Sampling["N_selected"] = row["N_models"]
                self.Sampling["N_sample_A"] = row["N_A"]
                self.Sampling["N_sample_B"] = row["N_B"]
                break

        # Save file
        columns = (
            ["Total_Score"]
            + [
                x for x in S_clusters.columns.values
                if x not in ["cluster", "Total_Score", "N_models",
                             "N_A", "N_B"]
            ]
            + ["N_models", "N_A", "N_B"]
        )
        S_clusters = S_clusters[columns]
        S_clusters.to_csv(
            os.path.join(self.analysis_dir, "summary_hdbscan_clustering.dat"),
            index=True,
        )
        print("Clustering summary: ")
        print(S_clusters)

    def plot_hdbscan_runs_info(self, S_comb):
        """
        Plot the number of models that come from each run
        """

        # Runs in each half
        runs_A = []
        for t in self.dir_halfA:
            runs_A.append([x for x in t.split("/") if self.dir_name in x][0])
        runs_A.sort(key=lambda x: float(x.split(self.dir_name)[-1]))
        runs_B = []
        for t in self.dir_halfB:
            runs_B.append([x for x in t.split("/") if self.dir_name in x][0])
        runs_B.sort(key=lambda x: float(x.split(self.dir_name)[-1]))

        # Count models per-run
        clusters = list(set(S_comb["cluster"]))
        clusters = [cl for cl in clusters if cl >= 0]
        if len(clusters) == 0:
            clusters = [-1]

        for cl in clusters:
            counts_traj_A = pd.DataFrame(columns=["traj", "counts_A"])
            counts_traj_B = pd.DataFrame(columns=["traj", "counts_B"])
            counts_traj_A["traj"] = runs_A
            counts_traj_B["traj"] = runs_B

            HH_cluster = S_comb[S_comb["cluster"] == cl]
            # Select two-halves
            HA = HH_cluster[(HH_cluster["half"] == "A")]
            HB = HH_cluster[(HH_cluster["half"] == "B")]
            if len(HA["half"].values) >= 0 and len(HB["half"].values) >= 0:
                for traj in counts_traj_A["traj"].values:
                    counts_traj_A.loc[counts_traj_A["traj"] == traj,
                                      "counts_A"] = len(HA[HA["traj"] == traj])
                for traj in counts_traj_B["traj"].values:
                    counts_traj_B.loc[counts_traj_B["traj"] == traj,
                                      "counts_B"] = len(HB[HB["traj"] == traj])
                n_max = np.max(
                    list(counts_traj_A["counts_A"].values)
                    + list(counts_traj_B["counts_B"].values)
                )

                # Plot
                fig, ax = pl.subplots(figsize=(10.0, 4.0), nrows=1, ncols=2)
                axes = ax.flatten()

                axes[0].bar(
                    np.arange(len(runs_A)),
                    counts_traj_A["counts_A"],
                    color="gold",
                    label=runs_A,
                )
                axes[0].set_title("Models sample A")
                axes[0].set_xlabel("Run number")
                axes[0].set_ylabel("Number of selected models")
                axes[0].set_xticks(np.arange(0, len(runs_A), 1))
                axes[0].set_xticklabels(runs_A, rotation=45)
                axes[0].set_ylim([0, n_max])

                axes[1].bar(
                    np.arange(len(runs_B)),
                    counts_traj_B["counts_B"],
                    color="orangered",
                    label=runs_B,
                )
                axes[1].set_title("Models sample B")
                axes[1].set_xlabel("Run number")
                axes[1].set_ylabel("Number of selected models")
                axes[1].set_xticks(np.arange(0, len(runs_B), 1))
                axes[1].set_xticklabels(runs_B, rotation=45)
                axes[1].set_ylim([0, n_max])

                pl.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.5)
                fig.savefig(
                    os.path.join(
                        self.analysis_dir,
                        "plot_run_models_cluster" + str(cl)
                        + "." + self.plot_fmt,
                    )
                )
                pl.close()

    def write_hdbscan_clustering(self, S_comb):

        """
        Write the frames information for each cluster
        """

        # Remove files from previous runs
        try:
            os.remove(os.path.join(self.analysis_dir,
                                   "selected_models_A_cluster*"))
            os.remove(os.path.join(self.analysis_dir,
                                   "selected_models_B_cluster*"))
        except:  # noqa: E722
            pass

        print("Selecting and writing models to extract ...")
        clusters = list(set(S_comb["cluster"]))
        clusters = [cl for cl in clusters if cl >= 0]

        if len(clusters) == 0:
            clusters = [-1]

        S_comb.loc[:, "frame_RMF3"] = S_comb.apply(
            lambda row: "h1_" + row.traj + "_"
            + str(int(row.MC_frame)) + ".rmf3"
            if row.half == "A"
            else "h2_" + row.traj + "_" + str(int(row.MC_frame)) + ".rmf3",
            axis=1,
        )

        clus_sel = 0
        for cl in clusters:
            HH_cluster = S_comb[S_comb["cluster"] == cl]

            # Select two-halves
            HA = HH_cluster[(HH_cluster["half"] == "A")]
            HB = HH_cluster[(HH_cluster["half"] == "B")]

            if len(HA["half"].values) >= 10 and len(HB["half"].values) >= 10:
                clus_sel += 1
                n = self.plot_scores_distributions(HA, HB, cl)

                # Write to csv
                HA.to_csv(
                    os.path.join(
                        self.analysis_dir,
                        "selected_models_A_cluster" + str(cl)
                        + "_detailed.csv",
                    )
                )
                HB.to_csv(
                    os.path.join(
                        self.analysis_dir,
                        "selected_models_B_cluster" + str(cl)
                        + "_detailed.csv",
                    )
                )

                # Select n model from
                if (
                    int(n) > self.number_models_out
                    or len(HH_cluster) > self.number_models_out
                ):
                    HH_sel = HH_cluster.sample(n=self.number_models_out)
                    HH_sel_A = HH_sel[(HH_sel["half"] == "A")]
                    HH_sel_B = HH_sel[(HH_sel["half"] == "B")]
                    HH_sel_A.to_csv(
                        os.path.join(
                            self.analysis_dir,
                            "selected_models_A_cluster"
                            + str(cl)
                            + "_detailed_random.csv",
                        )
                    )
                    HH_sel_B.to_csv(
                        os.path.join(
                            self.analysis_dir,
                            "selected_models_B_cluster"
                            + str(cl)
                            + "_detailed_random.csv",
                        )
                    )

        if clus_sel == 0:
            print(
                "WARNING: No models were selected because the "
                "simulations are not converged."
            )

    def plot_hdbscan_clustering(self, S_comb_sel, selected_scores):
        print("Generating HDBSCAN clustering plot ...")

        palette = color_palette[: len(np.unique(S_comb_sel["cluster"]))]
        cluster_colors = [palette[col] for col in S_comb_sel["cluster"]]

        n_sel = len(selected_scores)
        fig = pl.figure(figsize=(2 * n_sel, 2 * n_sel))
        gs = gridspec.GridSpec(n_sel, n_sel)
        i = 0
        for s1, s2 in itertools.product(selected_scores, repeat=2):
            ax = pl.subplot(gs[i])
            if s1 == s2:
                ax.hist(S_comb_sel[s1], 20, histtype="step",
                        color="b", alpha=0.5)
            else:
                ax.scatter(
                    S_comb_sel[s2],
                    S_comb_sel[s1],
                    c=np.array(cluster_colors),
                    s=3.0,
                    alpha=0.3,
                )

            i += 1

        # Add horizontal labels (top plots)
        for i in range(n_sel):
            fig.get_axes()[i].set_title(
                "%s" % (selected_scores[i]), fontsize=12)
        # Add vertical labels
        fig.get_axes()[0].set_ylabel("%s" % (selected_scores[0]), fontsize=12)
        k = 1
        for i in range(n_sel, n_sel * n_sel, n_sel):
            fig.get_axes()[i].set_ylabel(
                "%s" % (selected_scores[k]), fontsize=12)
            k += 1

        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(os.path.join(self.analysis_dir,
                                 "plot_clustering_scores.png"))
        pl.close()

    def do_extract_models(self, gsms_info, filename, gsms_dir):

        self.scores = self.manager.list()

        # Split the DF
        df_array = np.array_split(gsms_info, self.nproc)

        # Setup a list of processes that we want to run
        processes = [
            mp.Process(
                target=self.extract_models,
                args=(df_array[x], filename, gsms_dir)
            )
            for x in range(self.nproc)
        ]

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

        # Write scores to file
        np.savetxt(os.path.join(gsms_dir, filename + ".txt"),
                   np.array(self.scores))

    def write_GSMs_info(self, gsms_info, filename):
        gsms_info.to_csv(os.path.join(self.analysis_dir, filename),
                         index=False)

    def get_models_to_extract(self, f):
        # Get models to extract from file
        DD = pd.read_csv(f)
        return DD

    def get_sample_of_models_to_extract(self, file_A, file_B):
        # todo
        # DD_A = pd.read_csv(file_A)
        # DD_B = pd.read_csv(file_B)

        # Take samples and see when it plateaus
        pass

    def extract_models(self, gsms_info, filename, gsms_dir):
        """
        Use rmf_slice to extract the GSMs
        """
        for row in gsms_info.itertuples():
            ID = row.traj
            fr = int(row.MC_frame)
            fr_rmf = int(row.rmf_frame_index)
            rmf_file = row.rmf3_file
            traj_in = os.path.join(ID, rmf_file)
            file_out = os.path.join(
                gsms_dir, filename + "_" + str(ID) + "_" + str(fr) + ".rmf3"
            )

            os.system(
                "rmf_slice -q " + traj_in + " " + file_out
                + " --frame " + str(fr_rmf)
            )

            # Collect scores
            self.scores.append(row.Total_Score)

    def do_extract_models_single_rmf(
        self,
        gsms_info,
        out_rmf_name,
        traj_dir,
        analysis_dir,
        scores_prefix="scores",
        clean_rmfs=True,
        sel_state=0,
    ):
        """
        Extract all models to a single RMF file

        First, extract models from each independent run into an RMF file
        in that run directory (multithreaded)

        Then, concatenate these RMFs into one RMF in the analysis_dir
        along with the extracted scores from each file

        @param gsms_info : pandas dataframe with good scoring model information
        @param out_rmf_name : the filename for the output RMF
               (and "temporary" rmfs created in each trajectory directory)
        @param traj_dir : the folder containing all trajectory folders
        @param analysis_dir : the analysis folder for the output files
        @param scores_prefix : Prefix for the output scores file
        @param clean_rmfs : if True, delete the RMF files created in
               each trajectory directory
        """

        self.scores = self.manager.list()

        # Split the DF into pieces based on trajectory

        # Find the number of trajectories
        traj_dirs = gsms_info["traj"].unique()
        split_dfs = []  # Dataframe split into pieces
        filenames = []  # RMF file names
        scorefiles = []

        for td in traj_dirs:
            split_dfs.append(gsms_info[gsms_info["traj"] == td])
            filenames.append(os.path.join(traj_dir, td, out_rmf_name))
            scorefiles.append(os.path.join(traj_dir, td,
                                           scores_prefix + ".txt"))

        # Setup a list of processes that we want to run
        processes = [
            mp.Process(
                target=self.extract_models_to_single_rmf,
                args=(split_dfs[x], filenames[x], traj_dir,
                      scorefiles[x], sel_state),
            )
            for x in range(len(filenames))
        ]

        # Run processes
        for p in processes:
            p.start()

        # Exit the completed processes
        for p in processes:
            p.join()

        output_rmf = os.path.join(analysis_dir, out_rmf_name)
        output_score_file = os.path.join(analysis_dir, scores_prefix + ".txt")

        # check if rmf_cat is available
        if not shutil.which("rmf_cat"):
            raise RuntimeError("rmf_cat binary not found on path.")

        # concatenate RMF files
        print("Concatenating", len(filenames), "RMF files to", output_rmf)
        subprocess.check_call(["rmf_cat"] + filenames + [output_rmf])

        # concatenate score files
        with open(output_score_file, "wb") as outf:
            for scorefile in scorefiles:
                with open(scorefile, "rb") as scoref:
                    shutil.copyfileobj(scoref, outf)

        # Clean up our temporary files
        if clean_rmfs:
            for fn in filenames:
                os.remove(fn)
            for sf in scorefiles:
                os.remove(sf)

    def extract_models_to_single_rmf(
        self, inf, output_rmf, top_dir, scores_file, sel_state=0
    ):
        # Given a dataframe from get_models_to_extract(), extract these
        # models and place into a single RMF file: output_rmf

        scores = []
        i = 0

        # Initialize output RMF file
        row1 = inf.iloc[0]
        rmf_file = os.path.join(top_dir, row1.traj, row1.rmf3_file)
        # Create model and import hierarchies from one of the RMF files
        m = IMP.Model()
        f = RMF.open_rmf_file_read_only(rmf_file)
        h0 = IMP.rmf.create_hierarchies(f, m)[0]
        states = IMP.atom.get_by_type(h0, IMP.atom.STATE_TYPE)
        fh_out = RMF.create_rmf_file(output_rmf)
        for i, s in enumerate(states):
            if str(i) == sel_state:
                print("-------", str(i), sel_state)
                p = IMP.Particle(m, "System")
                hier_temp = IMP.atom.Hierarchy.setup_particle(p)
                hier_temp.add_child(s)
                IMP.rmf.add_hierarchy(fh_out, hier_temp)
        del f

        # Get a list of all the replica rmf3 files
        replicas = inf["rmf3_file"].unique()

        # Cycle through models by individual replica so we only open
        # one RMF at a time
        for rep in replicas:
            rep_inf = inf[inf["rmf3_file"] == rep]
            rmf_file = os.path.join(top_dir, row1.traj, rep)
            f = RMF.open_rmf_file_read_only(rmf_file)
            IMP.rmf.link_hierarchies(f, [h0])
            for (row_id, row) in rep_inf.iterrows():
                # t = row.traj
                # fr = int(row.MC_frame)
                fr_rmf = int(row.rmf_frame_index)

                # Sometimes individual frames don't print out for a number
                # of reasons. Just skip over these
                try:
                    IMP.rmf.load_frame(f, RMF.FrameID(fr_rmf))
                except:  # noqa: E722
                    continue
                IMP.rmf.save_frame(fh_out, str(i))

                if i % 1000 == 0:
                    print("Writing frame:", i, "of", len(inf.index),
                          "for", output_rmf)

                # Collect scores
                scores.append(row.Total_Score)
                i += 1
            del f

        # Write scores to file
        del fh_out
        np.savetxt(scores_file, np.array(scores))

    def create_gsms_dir(self, d):
        """
        Create directories for GSM.
        If already present, rename old one
        """
        if os.path.isdir(d):
            os.rename(d, "%s.old_%d" % (d, random.randint(0, 100)))
        os.mkdir(d)

    def plot_scores_distributions(self, HA, HB, cl):
        """
        Plot distribution of GSMs of both halves
        """

        n_bins = 20

        scores_A = HA["Total_Score"]
        scores_B = HB["Total_Score"]

        fig = pl.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(1, 2, width_ratios=[0.5, 0.5],
                               height_ratios=[1.0])

        min_score = np.min([scores_A.min(), scores_B.min()])
        max_score = np.max([scores_A.max(), scores_B.max()])
        # Plot distributions
        ax = pl.subplot(gs[0])
        ax.hist(
            scores_A,
            n_bins,
            histtype="step",
            stacked=True,
            fill=False,
            color="orangered",
        )
        ax.hist(
            scores_B, n_bins, histtype="step", stacked=True,
            fill=False, color="blue"
        )
        ax.set_xlabel("Total scores (a.u.)")
        ax.set_ylabel("Density")
        ax.set_title("Scores distribution")

        # Plot expected score from random set
        nn = len(scores_A) + len(scores_B)

        M = np.arange(
            int(min(len(scores_A), len(scores_B)) / 20.0),
            min(len(scores_A), len(scores_B)),
            int(nn / 20.0),
        )
        if len(M) <= 10.0:
            M = np.arange(
                int(min(len(scores_A), len(scores_B)) / 10.0),
                min(len(scores_A), len(scores_B)),
                int(min(len(scores_A), len(scores_B)) / 10.0),
            )

        RH1 = []
        RH2 = []
        for m in M[1:]:
            D1 = []
            D2 = []
            for d in range(20):
                D1.append(min(random.sample(list(scores_A), m)))
                D2.append(min(random.sample(list(scores_B), m)))
            RH1.append((m, np.mean(D1), np.std(D1)))
            RH2.append((m, np.mean(D2), np.std(D2)))

        RH1 = np.array(RH1)
        RH2 = np.array(RH2)

        hits = 0
        n = 0
        for s in range(len(RH1) - 1):
            dA = RH1[s, 1] - RH1[s + 1, 1]
            dB = RH2[s, 1] - RH2[s + 1, 1]
            if dA < RH1[s + 1, 2] and dB < RH2[s + 1, 1]:
                hits += 1
            if hits == 4:
                n = RH1[s + 1, 0]
                continue

        ax = pl.subplot(gs[1])
        ax.errorbar(RH1[:, 0], RH1[:, 1], yerr=RH1[:, 2],
                    c="orangered", fmt="o")
        ax.errorbar(RH2[:, 0], RH2[:, 1], yerr=RH2[:, 2], c="blue", fmt="o")
        ax.axvline(n, color="grey")
        ax.axhline(max_score, color="grey", ls="dashed", lw=3)
        ax.set_xlim([0.8 * M[0], 1.1 * M[-1]])
        ax.set_ylim([0.9 * min_score, 1.05 * max_score])

        ax.set_xlabel("Number of models")
        ax.set_ylabel("Minimum score (a.u.)")
        ax.set_title("Convergence of scores")

        pl.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.5)
        fig.savefig(
            os.path.join(
                self.analysis_dir,
                "plot_scores_convergence_cluster" + str(cl) + "."
                + self.plot_fmt,
            )
        )
        pl.close()

        return n

    def check_XLs_ambiguity(self, all_keys):
        """
        Input: DF with XLs distances
        Output: Dictionary of XLs that should be treated as ambiguous
        """
        xls_ids = {}

        if self.Multiple_XLs_restraints:
            for type_XLs in self.XLs_cutoffs.keys():
                xls_ids[type_XLs] = {}
                for i, xl in enumerate(all_keys):
                    if ("Distance_" in xl) and (type_XLs in xl):
                        vals = xl.split("|")
                        id = vals[2].split(".")[0]
                        # p1 = vals[3].split('.')[0]
                        # r1 = vals[4]
                        # p2 = vals[5].split('.')[0]
                        # r2 = vals[6]
                        if id in xls_ids[type_XLs].keys():
                            xls_ids[type_XLs][id].append(xl)
                        else:
                            xls_ids[type_XLs][id] = [xl]
        else:
            for i, xl in enumerate(all_keys):
                if "Distance_" in xl:
                    vals = xl.split("|")
                    id = vals[2].split(".")[0]
                    # p1 = vals[3].split('.')[0]
                    # r1 = vals[4]
                    # p2 = vals[5].split('.')[0]
                    # r2 = vals[6]
                    if id in xls_ids.keys():
                        xls_ids[id].append(xl)
                    else:
                        xls_ids[id] = [xl]
        return xls_ids

    def get_psi_stats(self):
        """
        Organize Psi values into DataFrame
        Get mean value for extracting models
        (Does not work for atomic XLs restraint)
        """
        key_0 = list(self.S_info_all.keys())[0]
        XLs_nuis = [
            k
            for k in self.S_info_all[key_0].columns.values
            if "CrossLinkingMassSpectrometryRestraint_Psi_" in k
        ]
        occams_nuis = [
            k
            for k in self.S_info_all[key_0].columns.values
            if ("OccamsRestraint" in k) and ("psi" in k)
        ]
        all_nuis = XLs_nuis + occams_nuis
        stats_nuis = sum([["mean_" + i, "std_" + i] for i in all_nuis], [])
        DF_stat_nuis = pd.DataFrame(columns=["traj"] + stats_nuis)

        for k, v in self.S_info_all.items():
            sel_nuis_mean = list(v[all_nuis].mean())
            sel_nuis_std = list(v[all_nuis].std())
            # DF_stat_nuis = DF_stat_nuis.append(
            #     pd.Series(
            #         [k] + sel_nuis_mean + sel_nuis_std,
            #         index=DF_stat_nuis.columns.values,
            #     ),
            #     ignore_index=True,
            # )
            DF_stat_nuis = pd.concat([DF_stat_nuis, pd.DataFrame(pd.Series(
                                [k] + sel_nuis_mean + sel_nuis_std,
                                index=DF_stat_nuis.columns.values).to_dict(),
                                index=[0])],
                            ignore_index=True)

        DF_stat_nuis.to_csv(os.path.join(self.analysis_dir,
                                         "Stat_all_nuisances.csv"))

    def analyze_trajectory_XLs(self, S_dist, S_info, atomic_XLs, traj_number,
                               ts_max):
        sel_XLs_nuis = [
            v
            for v in S_info.columns.values
            if "CrossLinkingMassSpectrometryRestraint_Psi" in v
        ]
        DF_XLs_psi = S_info[sel_XLs_nuis]
        psi_head = self.get_str_match(DF_XLs_psi.columns.values)
        # Get XLS satisfaction, append to S_dist
        XLs_satif_fields = []
        if self.Multiple_XLs_restraints:
            for type_XLs in self.XLs_cutoffs.keys():
                XLs_satif = self.get_XLs_satisfaction(
                    S_dist, atomic_XLs, type_XLs=type_XLs
                )
                temp_name = "XLs_satif_" + type_XLs.rstrip()
                S_info = S_info.assign(XLs_satif=pd.Series(XLs_satif))
                S_info.rename(columns={"XLs_satif": temp_name}, inplace=True)
                XLs_satif_fields.append(temp_name)
        elif self.Multiple_psi_values:
            all_psis = [
                v.split(psi_head)[1]
                for v in DF_XLs_psi.columns.values[1:]
                if "std" not in v
            ]
            for type_psi in all_psis:
                XLs_satif = self.get_XLs_satisfaction(
                    S_dist, atomic_XLs, type_psi=type_psi
                )
                temp_name = "XLs_satif_" + type_psi
                S_info = S_info.assign(XLs_satif=pd.Series(XLs_satif))
                S_info.rename(columns={"XLs_satif": temp_name}, inplace=True)
                XLs_satif_fields.append(temp_name)
        else:
            XLs_satif = self.get_XLs_satisfaction(S_dist, atomic_XLs)
            S_info = S_info.assign(XLs_satif=pd.Series(XLs_satif))
            XLs_satif_fields.append("XLs_satif")

        file_out_xls = "plot_XLs_%s.%s" % (traj_number, self.plot_fmt)
        self.plot_XLs_satisfaction(S_info, ts_max, file_out_xls)

        return S_info

    def get_XLs_satisfaction(self, S_dist, atomic_XLs, type_XLs=None,
                             type_psi=None):
        if type_XLs and not type_psi:
            dist_columns = [
                x for x in S_dist.columns.values
                if ("Distance_" in x and type_XLs in x)
            ]
            cutoff = self.XLs_cutoffs[type_XLs]
        elif type_psi and not type_XLs:
            dist_columns = [
                x for x in S_dist.columns.values
                if ("Distance_" in x and type_psi in x)
            ]
            cutoff = list(self.XLs_cutoffs.values())[0]
        elif type_XLs and type_psi:
            dist_columns = [
                x
                for x in S_dist.columns.values
                if ("Distance_" in x and type_XLs in x and type_psi in x)
            ]
            cutoff = self.XLs_cutoffs[type_XLs]
        elif atomic_XLs:
            dist_columns = [x for x in S_dist.columns.values
                            if ("BestDist" in x)]
            cutoff = list(self.XLs_cutoffs.values())[0]
        else:
            dist_columns = [x for x in S_dist.columns.values
                            if "Distance_" in x]
            cutoff = list(self.XLs_cutoffs.values())[0]

        # Only distance columns
        XLs_dists = S_dist[dist_columns]

        # Check for ambiguity
        if self.ambiguous_XLs_restraint is True:
            ambiguous_XLs_dict = self.check_XLs_ambiguity(
                S_dist.columns.values)
            min_XLs = pd.DataFrame()
            if self.Multiple_XLs_restraints:
                for k, v in ambiguous_XLs_dict[type_XLs].items():
                    min_XLs[k] = XLs_dists[v].min(axis=1)
            else:
                for k, v in ambiguous_XLs_dict.items():
                    min_XLs[k] = XLs_dists[v].min(axis=1)
            perc_per_step = list(
                min_XLs.apply(lambda x: sum(x <= cutoff) / len(x), axis=1)
            )

        else:
            perc_per_step = list(
                XLs_dists.apply(lambda x: sum(x <= cutoff) / len(x), axis=1)
            )

        return perc_per_step

    def summarize_XLs_info(
        self, Multiple_XLs_restraints=False, ambiguous_XLs_restraint=False
    ):
        # Remove files from previous runs
        try:
            os.remove(os.path.join(self.analysis_dir, "plot_XLs_distances_*"))
            os.remove(os.path.join(
                self.analysis_dir, "XLs_satisfaction_cluster_*"))
            os.remove(os.path.join(
                self.analysis_dir, "XLs_distances_cluster_*"))
        except:  # noqa: E722
            pass

        # If re-renning analysis
        if self.Multiple_XLs_restraints:
            pass
        else:
            self.Multiple_XLs_restraints = Multiple_XLs_restraints
        try:
            if self.ambiguous_XLs_restraint:
                pass
        except:  # noqa: E722
            self.ambiguous_XLs_restraint = ambiguous_XLs_restraint

        unique_clusters = np.sort(
            list(set(self.S_comb_dist_clustering["cluster"])))
        print("Summarize XLs, unique_clusters", unique_clusters)
        if self.Multiple_XLs_restraints:
            for type_XLs in self.XLs_cutoffs.keys():
                cutoff = self.XLs_cutoffs[type_XLs]
                for cl in unique_clusters:
                    # Boxplot XLs distances
                    self.boxplot_XLs_distances(
                        cluster=cl, type_XLs=type_XLs, cutoff=cutoff
                    )
                    # XLs satisfaction data
                    self.get_XLs_details(cluster=cl, type_XLs=type_XLs)
        else:
            cutoff = list(self.XLs_cutoffs.values())[0]
            for cl in unique_clusters:
                # Boxplot XLs distances
                self.boxplot_XLs_distances(cluster=cl, cutoff=cutoff)
                # XLs satisfaction data
                self.get_XLs_details(cluster=cl)

        # XLs satisfaction data for all models
        self.get_XLs_details(cluster="All")

    def get_XLs_details(self, cluster=0, type_XLs=None):
        """
        For GSM, determine for each XLs how often it is satisfied.
        """
        if type_XLs:
            dist_columns = [
                v
                for v in self.S_comb_dist_clustering.columns.values
                if "Distance" in v and type_XLs in v
            ]
            cutoff = [v for k, v in self.XLs_cutoffs.items()
                      if type_XLs in k][0]
        else:
            dist_columns = [
                v for v in self.S_comb_dist_clustering.columns.values
                if "Distance" in v
            ]
            cutoff = list(self.XLs_cutoffs.values())[0]

        if cluster != "All":
            dXLs_cluster = self.S_comb_dist_clustering.loc[
                self.S_comb_dist_clustering["cluster"] == cluster, dist_columns
            ]
        else:
            dXLs_cluster = self.S_comb_dist_clustering.loc[:, dist_columns]

        stats_XLs = pd.DataFrame()
        stats_XLs["mean"] = dXLs_cluster.mean()
        stats_XLs["std"] = dXLs_cluster.std()
        stats_XLs["min"] = dXLs_cluster.min()
        stats_XLs["max"] = dXLs_cluster.max()

        stats_XLs["perc_satif"] = dXLs_cluster.apply(
            lambda x: float(len(x[x < cutoff])) / float(len(x)), axis=0
        )

        if type_XLs:
            stats_XLs.to_csv(
                os.path.join(
                    self.analysis_dir,
                    "XLs_satisfaction_"
                    + type_XLs
                    + "_cluster_"
                    + str(cluster)
                    + ".csv",
                )
            )
            dXLs_cluster.to_csv(
                os.path.join(
                    self.analysis_dir,
                    "XLs_distances_" + type_XLs + "_cluster_"
                    + str(cluster) + ".csv",
                )
            )
        else:
            stats_XLs.to_csv(
                os.path.join(
                    self.analysis_dir,
                    "XLs_satisfaction_cluster_" + str(cluster) + ".csv",
                )
            )
            dXLs_cluster.to_csv(
                os.path.join(
                    self.analysis_dir, "XLs_distances_cluster_"
                    + str(cluster) + ".csv"
                )
            )

    def summarize_sampling_info(self):
        """
        Save information regarding number of models
        for output table
        """

        self.Sampling["Number_of_replicas"] = list(
            set(self.Sampling["Number_of_replicas"])
        )[0]
        self.Sampling["N_total"] = np.sum(self.Sampling["N_total"])
        self.Sampling["N_equilibrated"] = np.sum(
            self.Sampling["N_equilibrated"])
        # self.Sampling['Replica_exchange_temperature_range'] = '1.0-3.0'

        DS = pd.DataFrame()
        for k, v in self.Sampling.items():
            DS[k] = [v]

        DS.to_csv(
            os.path.join(self.analysis_dir,
                         "summary_sampling_information.csv"),
            index=False,
        )

    def summarize_fit_to_information(self):
        return 0

    def update_mmcif(self, mmcif_file):
        return 0

    def plot_XLs_satisfaction(self, S_info, ts_max, file_out):
        c = ["gold", "orange", "red", "blue", "green"]
        n_bins = 20

        XLs_percent = [k for k in S_info.columns.values if "XLs_satif" in k]
        XLs_nuis = [
            k
            for k in S_info.columns.values
            if "CrossLinkingMassSpectrometryRestraint_Psi_" in k
        ]

        XLs_percent.sort()
        XLs_nuis.sort()

        fig, ax = pl.subplots(figsize=(10.0, 4.0), nrows=1, ncols=3)
        axes = ax.flatten()
        for i, c in enumerate(XLs_percent):
            label = c
            axes[0].plot(S_info["MC_frame"][::20], S_info[c].loc[::20],
                         label=label)
        axes[0].set_title("XLs restraint satisfaction", fontsize=14)
        axes[0].set_xlabel("Step", fontsize=12)
        axes[0].set_ylabel("Percent Satisfied", fontsize=12)
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(handles[::-1], labels[::-1])

        for i, c in enumerate(XLs_nuis):
            label = c.split("CrossLinkingMassSpectrometryRestraint_")[-1]
            axes[1].plot(S_info["MC_frame"][::20], S_info[c].loc[::20],
                         label=label)
        axes[1].set_title("Psi nuisance parameters", fontsize=14)
        axes[1].set_xlabel("Step", fontsize=12)
        axes[1].set_ylabel("Psi", fontsize=12)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[1].legend(handles[::-1], labels[::-1])

        for i, c in enumerate(XLs_nuis):
            label = c.split("CrossLinkingMassSpectrometryRestraint_")[-1]
            axes[2].hist(
                S_info[c].loc[ts_max:], n_bins, histtype="step", fill=False,
                label=label)

        axes[2].set_title("Psi nuisance parameters", fontsize=14)
        axes[2].set_xlabel("Psi", fontsize=12)
        axes[2].set_ylabel("Density", fontsize=12)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[2].legend(handles[::-1], labels[::-1])

        pl.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.5)
        fig.savefig(os.path.join(self.analysis_dir, file_out))
        pl.close()

    def boxplot_XLs_distances(self, cluster=0, type_XLs=None, cutoff=30.0):
        if type_XLs:
            file_out = (
                "plot_XLs_distance_distributions_cl"
                + str(cluster)
                + "_"
                + str(type_XLs)
                + "."
                + self.plot_fmt
            )
        else:
            file_out = (
                "plot_XLs_distance_distributions_cl"
                + str(cluster)
                + "."
                + self.plot_fmt
            )

        if type_XLs:
            dist_columns = [
                x
                for x in self.S_comb_dist_clustering.columns.values
                if ("Distance_" in x) and (type_XLs in x)
            ]
        else:
            dist_columns = [
                x
                for x in self.S_comb_dist_clustering.columns.values
                if "Distance_" in x
            ]
        dXLs_cluster = self.S_comb_dist_clustering.loc[
            self.S_comb_dist_clustering["cluster"] == cluster, dist_columns
        ]

        dXLs_unique = pd.DataFrame()
        if self.ambiguous_XLs_restraint is True:
            ambiguous_XLs_dict = self.check_XLs_ambiguity(dist_columns)
            if type_XLs:
                for k, v in ambiguous_XLs_dict.items():
                    XLs_sele = dXLs_cluster.loc[:, v].mean()
                    XLs_min = XLs_sele.idxmin()
                    dXLs_unique[XLs_min] = dXLs_cluster[XLs_min]
            else:
                for k, v in ambiguous_XLs_dict.items():
                    XLs_sele = dXLs_cluster.loc[:, v].mean()
                    XLs_min = XLs_sele.idxmin()
                    dXLs_unique[XLs_min] = dXLs_cluster[XLs_min]
        else:
            dXLs_unique = dXLs_cluster

        # Get distances and order based on the mean
        columns = ["entry", "mean", "id"]
        XLs_ids = pd.DataFrame(columns=columns)
        min_all = []
        for i, v in enumerate(dXLs_unique.columns.values):
            m = np.mean(dXLs_unique[v])
            ll = v.split("_|")[1]
            label = "|".join(ll.split("|")[2:6])
            if m > 0:
                XLs_ids = pd.concat([XLs_ids, pd.DataFrame(
                                pd.Series([int(i), m, label],
                                    index=columns).to_dict(),
                                index=[0])],
                            ignore_index=True)
                min_all.append(m)
        XLs_ids = XLs_ids.sort_values(by=["mean"])
        labels_ordered = XLs_ids["id"].values

        # For plot layout
        S_sorted = np.array(dXLs_unique)[
            :, XLs_ids["entry"].values.astype("int")]
        n_xls = len(labels_ordered)
        n_plots = int(math.ceil(n_xls / 50.0))
        n_frac = int(math.ceil(n_xls / float(n_plots)))

        # Generate plot
        fig, ax = pl.subplots(figsize=(12, 6.0 * n_plots), nrows=n_plots,
                              ncols=1)
        if n_plots == 1:
            ax = [ax]
        for i in range(n_plots):
            if i == n_plots - 1:
                _ = ax[i].boxplot(
                    S_sorted[:, (i * n_frac): -1], patch_artist=True,
                    showfliers=False)
                ax[i].set_xticklabels(
                    labels_ordered[(i * n_frac): -1], rotation="vertical",
                    fontsize=10)
                max_y = np.max(S_sorted[:, -1]) + 25
                ax[i].set_ylim([0, (max_y - (max_y % 25))])
            else:
                _ = ax[i].boxplot(
                    S_sorted[:, (i * n_frac): ((i + 1) * n_frac)],
                    patch_artist=True,
                    showfliers=False,
                )
                ax[i].set_xticklabels(
                    labels_ordered[(i * n_frac): ((i + 1) * n_frac)],
                    rotation="vertical",
                    fontsize=10,
                )
                ax[i].set_ylim(
                    [0, np.max(S_sorted[:, (i + 1) * n_frac]) + 20.0])
            ax[i].axhline(y=cutoff, color="forestgreen", linestyle="-",
                          linewidth=1.5)
            ax[i].set_xticks(range(1, n_frac + 1))

            ax[i].set_ylabel("XLs distances (A)")
            ax[i].set_title("XLs distance distributions")

        pl.tight_layout()
        fig.savefig(os.path.join(self.analysis_dir, file_out))
        pl.close()

        # Plot histogram of best cluster distances
        if type_XLs:
            file_out_hist = (
                "plot_XLs_distance_histogram_cl"
                + str(cluster)
                + "_"
                + str(type_XLs)
                + "."
                + self.plot_fmt
            )
        else:
            file_out_hist = (
                "plot_XLs_distance_histogram_cl" + str(cluster) + "."
                + self.plot_fmt)
        self.plot_XLs_satisfaction_histogram(min_all, cutoff, file_out_hist)

    def plot_XLs_satisfaction_histogram(self, min_all, cutoff, file_out_hist):
        fig, ax = pl.subplots(figsize=(5, 5), nrows=1, ncols=1)

        ax.hist(min_all, 20, color="b", alpha=0.5)
        ax.axvline(x=cutoff, color="orange", alpha=0.7, lw=3)
        ax.set_xlabel("Distance (A)", fontsize=12)
        ax.set_ylabel("Number of XLs")
        ax.set_title("XLs satisfaction")
        pl.tight_layout()
        fig.savefig(os.path.join(self.analysis_dir, file_out_hist))
        pl.close()

    def plot_pEMAP_satisfaction(self, S_info, file_out):
        fig, ax = pl.subplots(figsize=(4.0, 4.0), nrows=1, ncols=1)

        pemap_satif = [
            v for v in S_info.columns.values
            if "pEMapRestraint_satisfaction" in v][0]

        ax.plot(
            S_info["MC_frame"].iloc[::10],
            S_info[pemap_satif].iloc[::10],
            color="orangered",
            alpha=0.8,
        )
        ax.set_title("pE-MAP restraint", fontsize=14)
        ax.set_xlabel("Step", fontsize=12)
        ax.set_ylabel("Percent satisfied", fontsize=12)

        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(os.path.join(self.analysis_dir, file_out))
        pl.close()

    def plot_Occams_satisfaction(self, Occams_info, file_out):
        """Plot percent of restraint satisfied and the distribution of
        the nuisances"""

        c = ["gold", "red", "blue", "green"]

        fig, ax = pl.subplots(figsize=(12.0, 4.0), nrows=1, ncols=3)
        # Occams satisfaction
        occams_satif = [
            k for k in Occams_info.columns.values
            if "OccamsRestraint_satisfied" in k]

        # Occams Psi nuisances
        occams_psi_nuis = [
            k
            for k in Occams_info.columns.values
            if ("OccamsRestraint" in k) and ("psi" in k)
        ]

        # Occams Sigma nuisances
        occams_sigma_nuis = [
            k
            for k in Occams_info.columns.values
            if ("OccamsRestraint" in k) and ("sigma" in k)
        ]

        for i, k in enumerate(occams_satif):
            label = k
            ax[0].plot(
                Occams_info["MC_frame"].loc[::20],
                Occams_info[k].loc[::20],
                alpha=0.8,
                label=label,
            )
        ax[0].set_title("Occams restraint", fontsize=14)
        ax[0].set_xlabel("Step", fontsize=12)
        ax[0].set_ylabel("Percent satisfied", fontsize=12)
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(handles[::-1], labels[::-1])

        for i, k in enumerate(occams_psi_nuis):
            label = k
            ax[1].plot(
                Occams_info["MC_frame"].loc[::20],
                Occams_info[k].loc[::20],
                alpha=0.8,
                label=label,
            )
        ax[1].set_title("Occams restraint nuisances", fontsize=14)
        ax[1].set_xlabel("Step", fontsize=12)
        ax[1].set_ylabel("Psi value", fontsize=12)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[1].legend(handles[::-1], labels[::-1])

        for i, k in enumerate(occams_sigma_nuis):
            label = k
            ax[2].plot(
                Occams_info["MC_frame"].loc[::10],
                Occams_info[k].loc[::10],
                color=c[i],
                alpha=0.8,
                label=label,
            )
        ax[2].set_title("Occams restraint nuisances", fontsize=14)
        ax[2].set_xlabel("Step", fontsize=12)
        ax[2].set_ylabel("Sigma value", fontsize=12)
        handles, labels = ax[2].get_legend_handles_labels()
        ax[2].legend(handles[::-1], labels[::-1])

        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(os.path.join(self.analysis_dir, file_out))
        pl.close()

    def substrings(self, s):
        for i in range(len(s)):
            for j in range(i, len(s)):
                yield s[i: j + 1]

    def get_str_match(self, strs):
        if len(strs) > 1:
            intersect = (set(self.substrings(strs[0]))
                         & set(self.substrings(strs[1])))
            return max(intersect, key=len)
        else:
            return strs[0]
