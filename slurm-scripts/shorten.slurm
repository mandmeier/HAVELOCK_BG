#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=16gb
#SBATCH --job-name=T145
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/shorten_5_out_traits_150_lowN_201104-185415/T145_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/shorten_5_out_traits_150_lowN_201104-185415/T145_%J.out

module load python/3.8

python /work/jyanglab/mmeier/HAVELOCK_BG/profiling/8.GWAS/shorten.py /work/jyanglab/mmeier/HAVELOCK_BG/largedata/GWAS/out_traits_150_lowN_201104-185415/T145.assoc.txt 1e-05 /work/jyanglab/mmeier/HAVELOCK_BG/largedata/GWAS/short_5_out_traits_150_lowN_201104-185415/T145.assoc.txt