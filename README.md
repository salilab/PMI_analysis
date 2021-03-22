These scripts are used to analyze a series of trajectories, look at their scores, determine when each score has equilibrated, and select models for further analysis. Sample scrips are in the *example* folder. These scripts assume that you are a series of IMP runs named run_0, run_1, ...

## Dependencies:
* multiprocessing
* hdbscan
* numpy
* pandas
* matplotlib
* pyrmsd

## Analysis steps:

1. To initialize the analysis class:

```python
AT = AnalysisTrajectories(out_dirs,
                          dir_name = 'run_',
                          analysis_dir = analysis_dir,
                          detect_equilibration=True,
                          burn_in_fraction=0.02,
                          nskip=20,
                          nproc=1,
                          plot_fmt='pdf')
```
The minimum set of arguments to specify are ```dir_name```: the prefix of the run folders, ```analysis_dir```: the directory where all analysis output will be written, and ```nproc```: the number of processors that are used to do parallel I/O. Additional options are: ```detect_equilibration``` which when set to ```False``` turns off rigorous statistical determination of the burn-in  or equilibration phase and instead just discards a ```burn_in_fraction``` (default 2% of the trajectory size) fraction of frames from the beginning of each independent run, to calculate statistics. ```nskip``` is used to *thin* out the samples collected for analysis, i.e. every ```nskip```^th  frame is selected. Finally, ```plot_fmt``` specifies the file extension of all matplotlib figures generated. You should only set ```detect_equilibration=False``` if you are testing your analysis script and need it to run fast without waiting for calculating the exact length of the monte carlo burn-in phase. 

2. Add the restraints that you want to be analyzed:

```python
XLs_cutoffs = {'DSSO':30.0}
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
```

If you set up more than one XLs restraint, you should include the *Multiple_XLs_restraints = True* flag in the analysis and include a cutoff for all of them in the *XLs_cutoffs* dictionary. For example, if you divided the DSSO XLs into a intra-subunit and inter-subunit datasets, there should be two elements in the *XLs_cutoffs* dictionary (even if they have the same cutoff):

```python
XLs_cutoffs = {'DSSO_Inter':30.0, 'DSSO_Intra':30.0}
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
                             Multiple_XLs_restraints = True)
```
Using the labels you used to setup the restraint and keys.

Similarly, if there is ambiguity in the XLs assignments (i.e. you have multiple copies of the same protein), you should use the *ambiguous_XLs_restraint = True* option:

```python
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
                             ambiguous_XLs_restraint = True)
```

If you want to analyze a restraint that is not standard to IMP, you can do so by including the handle associated to that restraint in the stat file. E.g. if your restraint is called ```COMDistanceRestraint_data_Score``` in the stat files produced after sampling (i.e. the restraint *handle*) and you'd like the restraint score to be referred to by a short name, say ```COM``` throughout the various analysis output data-frames:

```python
AT.set_analyze_score_only_restraint(handle='COMDistanceRestraint_data_Score',
                                    short_name='COM',
                                    do_sum=True)
```

The ```do_sum=True``` option will add up all restraint scores that start with the restraint handle as specified (e.g. ```COMDistanceRestraint_data_Score_1```, ```COMDistanceRestraint_data_Score_2```, etc). In general it may be desirable to add scores of the same type prior to score based clustering, since that decreases the dimensionality of the score space and helps the clustering process converge better. E.g., restraints like connectivity and excluded volume have their scores added up by default. Note that this method will only analyze the scores for the restraint, not nuisances or other associated statistics. 


3. Read the stat files to obtain the scores, nuisances parameters, and information about the rmf3 files:

```
AT.read_stat_files()
```

In this step we also automatically determine the equilibration time for each restraint and nuisance particles. Only models after equilibration are considered for further analysis (unless, as mentioned before the ```detect_equilibration``` option is set to False. )

Reference for equilibration detection:

Utilities for automatically detecting equilibrated region of molecular simulations. DOI: 10.1021/acs.jctc.5b00784

John D. Chodera <john.chodera@choderalab.org>

4. Obtain the statistics of the XLs restraint Psi nuisance parameter use:
```
AT.get_psi_stats()
```

5. Do HDBSCAN clustering for selected scores and/or nuisance parameters:
```
AT.hdbscan_clustering(['EV_sum', 'XLs_sum'])
```

If you added a non-standard restraint, you can use it in clustering with the name 'COMDistanceRestraint' previously used (everything before the first _)

6. Get information about XLs satisfaction:
```
AT.summarize_XLs_info()

```
This will create a series of files and plots summarizing the XLs distances and satisfaction in all the clusters obtained in step 6. Files ```plot_run_models_cluster*.<plot_fmt>``` show the number of models from each run that are in each cluster. Files ```plot_scores_convegence_cluster*.<plot_fmt>``` show the scores distribution for each cluster.

Files ```all_info_*.csv``` contain the information about all models after equilibration. These files can be used to re-run the clustering step:

```python
AT.read_models_info()
AT.hdbscan_clustering(['EV_sum', 'XLs_sum'])
```

After clustering a series of files are written with the information of frames in each cluster.

7. To re-rerun the clustering step without having to read all the stat files again, you can read the relevant information from the `all_info_*.csv` files:

```
XLs_cutoffs = {'DSSO_Inter':30.0, 'DSSO_Intra':30.0}
AT = AnalysisTrajectories(out_dirs,
                          analysis_dir = analys_dir,
                          nproc=nproc)
			 
AT.read_models_info(XLs_cutoffs)
AT.hdbscan_clustering(['EV_sum', 'XLs_sum'],
                        min_cluster_size=200,
                        min_samples=5,
                        skip=5)
AT.summarize_XLs_info()
			
```

Reference for HDBSCAN clustering:

http://hdbscan.readthedocs.io/en/latest/index.html

8. To extract the models from the rmf3 file, use the script `run_extract_models.py`. 

9. To test for convergence and do structural clustering, use the script `run_clustering.sh`

Reference for convergence and clustering anlysis:

Assessing exhaustiveness of stochastic sampling for integrative modeling of macromolecular structures. https://doi.org/10.1016/j.bpj.2017.10.005


10. If you know the structure of complex (i.e. you are benckmarking a method), you can determine the accuracy of the structural models using `get_accuracy_rmfs.py`
