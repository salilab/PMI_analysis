import IMP
import IMP.pmi
import IMP.pmi.macros
import RMF
import numpy as np


def get_coordinates_alignment(hier, selection=None):
    """Get coordinates from RMF frame"""

    coord_dict = {}

    if selection:
        for k, v in selection.items():
            sel = IMP.atom.Selection(
                hier,
                molecule=v[0],
                residue_indexes=np.arange(v[1], v[2], 1),
                resolution=IMP.atom.ALL_RESOLUTIONS,
                copy_index=v[3],
            ).get_selected_particles()
            coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in sel]
            coord_dict[k] = coords

    else:
        mols = IMP.pmi.tools.get_molecules(hier)
        print(mols)
        for m in mols:
            sel = IMP.atom.Selection(
                hier,
                molecule=m.get_name(),
                copy_index=IMP.atom.Copy(m).get_copy_index(),
                resolution=IMP.atom.ALL_RESOLUTIONS,
            ).get_selected_particles()

            coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in sel]

            coord_dict[m.get_name()] = coords

    return coord_dict


def transform_coordinates(hier, transformation):
    """Transform all coordinates"""
    rbs, beads = IMP.pmi.tools.get_rbs_and_beads(hier)
    for rb in rbs:
        IMP.core.transform(rb, transformation)
    for p in beads:
        temp_coord = IMP.core.XYZ(p)
        IMP.core.transform(temp_coord, transformation)


def get_reference_coordinates(rmf_in, selection=None):
    """Get reference coordiantes from frame 0"""
    m = IMP.Model()

    f = RMF.open_rmf_file_read_only(rmf_in)
    hier = IMP.rmf.create_hierarchies(f, m)[0]

    IMP.rmf.load_frame(f, RMF.FrameID(0))

    # Get coordinates from frame 0
    ref_coord = get_coordinates_alignment(hier, selection)
    del m, f
    return ref_coord


################################
def align_rmf(rmf_in, rmf_out,
              ref_coord, selection=None,
              frames=None, sel_state=0):

    fh_out = RMF.create_rmf_file(rmf_out)

    m = IMP.Model()
    f = RMF.open_rmf_file_read_only(rmf_in)
    print("Number of frames", f.get_number_of_frames())

    if not frames:
        frames = np.arange(0, f.get_number_of_frames(), 100)

    hier = IMP.rmf.create_hierarchies(f, m)[0]
    states = IMP.atom.get_by_type(hier, IMP.atom.STATE_TYPE)
    for i, s in enumerate(states):
        if i == sel_state:
            p = IMP.Particle(m, "System")
            hier_temp = IMP.atom.Hierarchy.setup_particle(p)
            hier_temp.add_child(s)
            IMP.rmf.add_hierarchy(fh_out, hier_temp)

    RMSD = []
    for i in frames:
        if i % 100 == 0:
            print("Frame:", i)
        IMP.rmf.load_frame(f, RMF.FrameID(i))

        temp_coord = get_coordinates_alignment(hier, selection)

        ali = IMP.pmi.analysis.Alignment(ref_coord, temp_coord)
        (rmsd, transformation) = ali.align()
        RMSD.append(rmsd)

        transform_coordinates(hier, transformation)
        IMP.rmf.save_frame(fh_out, str(i))

        del temp_coord

    del f

    # Save all RMSD values after alignment
    RMSD = np.array(RMSD)
    if "/" in rmf_in:
        name_in = rmf_in.split("/")[-1].split(".")[0]
    else:
        name_in = rmf_in.split(".")[0]
    out_RMSD = f"{RMSD}_{name_in}.txt"
    np.savetxt(out_RMSD, RMSD)
    print("Mean RMSD:", np.mean(RMSD))
