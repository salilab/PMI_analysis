#!/bin/bash
module load imp
module load python3/pyrmsd

export cl=1
export st=0

cp ../A_models_clust${cl}_${st}.txt ../scoresA.txt
cp ../B_models_clust${cl}_${st}.txt ../scoresB.txt
cp ../../density_DDI1_DDI2.txt density.txt
nohup imp_sampcon exhaust --sysname test --p ../ \
                          --rmfA A_models_clust${cl}_${st}_aligned.rmf3 \
                          --rmfB B_models_clust${cl}_${st}_aligned.rmf3 \
                          --density density.txt --m cpu_omp --cores 8 \
                          --gridsize 2.0 --voxel 3.0  --density_threshold 10 \
                          > clustering.log &


