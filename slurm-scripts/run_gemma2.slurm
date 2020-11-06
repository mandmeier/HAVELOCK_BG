#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --job-name=gemma_150_lowN
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/gemma_150_lowN_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/gemma_150_lowN_%J.out

module load python/3.8

python profiling/8.GWAS/run_gemma.py traits_150_lowN.txt
