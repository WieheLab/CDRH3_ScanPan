#!/bin/bash
#SBATCH --job-name=count
#SBATCH --out count-%j.out
#SBATCH --ntasks=1
#SBATCH --mem=35G
#SBATCH --partition=dhvi
#SBATCH --ntasks-per-node=1

for i in $(seq 2 22)
do
    echo "column position" ${i}

    cut -f ${i} igor_dist.txt | sort | uniq -c > counts_${i}.txt
done
