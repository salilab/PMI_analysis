from contact_maps import CMTable


#####################################################
# calculate contact frequency
#####################################################
CM = CMTable(out_dir='CMs_cluster0', GSMs_dir='a', clustering_dir='../RMF/',
             cluster=0, number_of_models=3,    cutoff=20.0, nproc=20)
CM.compute_contact_maps()

CM.get_close_contacts()
CM.plot_contact_maps_subunits()
