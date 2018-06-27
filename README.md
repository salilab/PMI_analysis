These scripts are used to analyze a series of trajactories, look at their scores, determine when each score has equilibrated, and select models for further analysis. Sample scrips are in the *example* folder. These scripts assume that you are a series of IMP runs named run_0, run_1, ...


1. To initialize the analysis class:

```
AT = AnalysisTrajectories(out_dirs,
                          analysis_dir = analys_dir,
                          nproc=nproc)

```

2. Add the restraints that you want to be analyzed:

```
XLs_cutoffs = {'DSSO':30.0}
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
```

3. Identify where each of the restraints scores and other relevant information is stored in the stat files

```
AT.get_restraint_fields()
```

4. Read the stat files to obtain the scores, nuisances parameters, and information about the rmf3 files:

```
AT.read_stat_files()
```

5. Obtain the statistics of the XLs restraint Psi nuisance parameter use:
```
AT.get_Psi_stats()
```

6. Do HDBSCAN clustering for selected scores and/or nuisance parameters:
```
AT.do_hdbscan_clustering(['EV_sum', 'XLs_sum', 'Psi_vals_0.01', 'Psi_vals_0.1'])
```

7. Get information about XLs satisfaction:
```
AT.get_XLs_details()

```

Files all_info_*.csv contain the information about all models after equilibration. These files can be used to re-run the clustering step:

```
AT.read_models_info()
AT.do_hdbscan_clustering(['EV_sum', 'XLs_sum'])
```

After clustering a series of files are written with the information of frames in each cluster.

To extract the models from the rmf3 file, use the script `run_extract_models.py`.