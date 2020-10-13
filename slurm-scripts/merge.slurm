#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=1500gb
#SBATCH --job-name=merge
#SBATCH --error=./LOG/merge_%J.err
#SBATCH --output=./LOG/merge_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.8_merge.R
