module load imp
module load python3/pyrmsd

export cl=1
export st=1

cp ../A_models_clust${cl}_${st}.txt ../scoresA.txt
cp ../B_models_clust${cl}_${st}.txt ../scoresB.txt
cp ../../density.txt density.txt
nohup imp-sampcon exhaust \
      -n test -p ../ -ra A_models_clust${cl}_${st}_aligned.rmf3 -rb B_models_clust${cl}_${st}_aligned.rmf3 -d density.txt \
      -m cpu_omp -c 8 -g 2.0 --voxel 3.0 --density_threshold 10 --prism > clustering.log &


