from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import RMF
import glob
import numpy as np
import pandas as pd
import os
import multiprocessing

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pylab as pl  # noqa: E402

mpl.rcParams.update({"font.size": 8})


def worker(task):
    method, args = task
    return method(*args)


class ParallelProcessor:
    def __init__(self):
        self.manager = multiprocessing.Manager()

    def parallel_process(self, tasks):
        # Use a Pool for multiprocessing
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(worker, tasks)
        return results


class AccuracyModels(object):
    def __init__(
        self,
        selection_dictionary,
        clustering_dir,
        rmf3_ref,
        rmf3_A,
        rmf3_B,
        out_header="all",
        skip=1,
    ):
        self.selection_dictionary = selection_dictionary
        self.clustering_dir = clustering_dir
        self.rmf3_ref = rmf3_ref
        self.rmf3_A = rmf3_A
        self.rmf3_B = rmf3_B
        self.out_header = out_header
        self.skip = int(skip)

        if "selection" not in self.selection_dictionary.keys():
            raise ValueError(
                'selection_dictionary need to include a selection key: selection_dictionary = {"selection"=["p1","p2"]} or {"selection"=[("r1","r2","p1")]}'
            )

        rhA = RMF.open_rmf_file_read_only(self.rmf3_A)
        self.nframes_A = rhA.get_number_of_frames()

        # Read clustering information
        self.read_cluster_frames()

        # Compute accuracies
        self.get_accuracy_clusters()

        # Plot and write results
        self.plot_accuracy_histograms()
        self.write_accuracy_values()

    def read_cluster_frames(self):
        # Dictionary cluster: frames
        self.frames_cluster = dict()

        # Read cluster files
        clusters = glob.glob(os.path.join(self.clustering_dir, "cluster.*.all.txt"))

        for cl in np.arange(0, len(clusters)):
            print("Reading cluster:", cl)
            frames_A = []
            frames_B = []
            for line in open(
                os.path.join(self.clustering_dir, f"cluster.{cl}.sample_A.txt"), "r"
            ):
                frames_A.append(int(line.split("\n")[0]))
            for line in open(
                os.path.join(self.clustering_dir, f"cluster.{cl}.sample_B.txt"), "r"
            ):
                frames_B.append(int(line.split("\n")[0]) - self.nframes_A)
            self.frames_cluster[cl] = {"A": frames_A, "B": frames_B}
            print(
                "Number of frames per cluster (A,B):",
                cl,
                len(self.frames_cluster[cl]["A"]),
                len(self.frames_cluster[cl]["B"]),
            )

    def get_accuracy_clusters(self):
        processor = ParallelProcessor()
        manager = processor.manager

        self.accuracy_values = manager.dict()
        for cluster in self.frames_cluster.keys():
            self.accuracy_values[cluster] = []

        for cluster, frames in self.frames_cluster.items():
            rmfs = [self.rmf3_A, self.rmf3_B]
            all_frames = [frames["A"], frames["B"]]

            tasks = [
                (
                    self.accuracy_rmfs,
                    (
                        cl,
                        rmf,
                        fr,
                    ),
                )
                for cl, rmf, fr in zip([cluster] * 2, rmfs, all_frames)
            ]

            processor.parallel_process(tasks)

    def accuracy_rmfs(self, cluster, rmf, frames):
        """
        Compute accuracy between set of frames and
        a reference structure
        """

        if self.skip > 1:
            frames = frames[:: self.skip]
        if len(frames) > 0:
            model = IMP.Model()
            pr = IMP.pmi.analysis.Precision(
                model, resolution=1, selection_dictionary=self.selection_dictionary
            )
            pr.set_precision_style("pairwise_rmsd")
            pr.add_structures(zip([rmf] * len(frames), frames), "set0")
            pr.set_reference_structure(self.rmf3_ref, 0)
            vals = pr.get_rmsd_wrt_reference_structure_with_alignment("set0", ("selection"))
            self.accuracy_values[cluster] += vals["selection"]["all_distances"]

    def write_accuracy_values(self):
        """
        Write accuracy values as well as
        the mean, min, and max
        """
        out_summary = open(
            os.path.join(
                self.clustering_dir, "accuracy_" + self.out_header + "_clusters.dat"
            ),
            "w",
        )

        out_summary.write("cluster mean min max N \n")
        for k, v in self.accuracy_values.items():
            out = open(
                os.path.join(
                    self.clustering_dir,
                    "accuracy_" + self.out_header + "_cl" + str(k) + ".dat",
                ),
                "w",
            )

            if len(np.array(v)) > 0:
                vv = np.array(v).astype(float)
                out_summary.write(
                    str(k)
                    + "\t"
                    + str(np.mean(vv))
                    + "\t"
                    + str(np.min(vv))
                    + "\t"
                    + str(np.max(vv))
                    + "\t"
                    + str(len(vv))
                    + "\n"
                )
                out.close()
        out_summary.close()

    def plot_accuracy_histograms(self):
        """
        Plot accuracy distribution for all clusters.
        Only cluster0 is shown in orange
        """

        colors = [
            "gold",
            "dodgerblue",
            "forestgreen",
            "palegreen",
            "darkviolet",
            "darkblue",
        ]
        n_bins = 10
        clusters_all = sorted(self.accuracy_values.keys())
        fig, ax = pl.subplots(figsize=(5.0, 5.0), nrows=1, ncols=1)
        for i, cluster in enumerate(clusters_all):
            if len(np.array(self.accuracy_values[cluster])) > 0:
                A = np.array(self.accuracy_values[cluster]).astype(float)
                if i == 0:
                    ax.hist(
                        A,
                        n_bins,
                        histtype="step",
                        fill=False,
                        color="orangered",
                        alpha=0.9,
                        lw=4,
                        label=f"cluster {cluster}",
                    )
                    ax.axvline(np.mean(A), color="orangered", alpha=0.9)
                else:
                    if len(clusters_all) <= 6:
                        color = colors[i - 1]
                    else:
                        color = "grey"
                    ax.hist(
                        A,
                        n_bins,
                        histtype="step",
                        fill=False,
                        color=color,
                        alpha=0.5,
                        label=f"cluster {cluster}",
                    )
                    ax.axvline(np.mean(A), color="grey", alpha=0.2, ls="dashed")

        if len(clusters_all) <= 6:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels)
        ax.set_title("Accuracy of selected models", fontsize=14)
        ax.set_xlabel("Accuracy (A)", fontsize=12)
        ax.set_ylabel("Number of models", fontsize=12)

        pl.tight_layout(pad=1.2, w_pad=1.5, h_pad=2.5)
        fig.savefig(
            os.path.join(
                self.clustering_dir, f"plot_accuracy_{self.out_header}_clusters.pdf"
            )
        )
