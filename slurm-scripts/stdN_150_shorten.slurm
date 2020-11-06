#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=16gb
#SBATCH --job-name=stdN_shorten
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/stdN_shorten_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/stdN_shorten_%J.out

module load python/3.8

python /work/jyanglab/mmeier/HAVELOCK_BG/profiling/8.GWAS/run_shorten.py out_traits_150_stdN_201104-164927 5


