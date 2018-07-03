module load cuda/9.0.176
module load imp-fast/last_ok_build

export top_dir=/home/ignacia/Research/pE-MAP/mod_xls_v2/
export analys_dir=$top_dir/analys
export mod_dir=$analys_dir/GSMs_cl0/
export name=RP_xls

cp $top_dir/density.txt $mod_dir
cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

ls -lta $analys_dir/GSMs_cl0/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster0_random.dat
ls -lta $analys_dir/GSMs_cl0/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster0_random.dat

nohup python ~/SOFTW/imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --mode cuda --cores 4 --align --density density.txt --gridsize 1.0 --gnuplot  --rmfs_A $analys_dir/selected_models_A_cluster0_random.dat --rmfs_B $analys_dir/selected_models_B_cluster0_random.dat > clustering.log
