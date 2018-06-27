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

Files all_info_*.csv contain the information about all models after equilibration. This files can be used to re-run the clustering step:

`AT.read_models_info()`
`AT.do_hdbscan_clustering(['EV_sum', 'XLs_sum'])`