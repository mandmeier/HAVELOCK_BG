#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=50gb
#SBATCH --job-name=snps
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/snps_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/snps_%J.out

ml python/3.8

python profiling/9.LD_haplotypes/subset_hapmap.py
