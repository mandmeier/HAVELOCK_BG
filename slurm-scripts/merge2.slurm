#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=500gb
#SBATCH --job-name=merge2
#SBATCH --error=./LOG/merge2_%J.err
#SBATCH --output=./LOG/merge2_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.8.2_merge2.R
