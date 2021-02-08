#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=2:00:00
#SBATCH --mem=50gb
#SBATCH --job-name=gwasmerge
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/gwasmerge_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/gwasmerge_%J.out

module load python/3.8

python /work/jyanglab/mmeier/HAVELOCK_BG/profiling/8.GWAS/merge.py short_1_out_traits_150_stdN_201109-103909
python /work/jyanglab/mmeier/HAVELOCK_BG/profiling/8.GWAS/merge.py short_1_out_traits_150_lowN_201109-100251
