#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=500gb
#SBATCH --job-name=merge3
#SBATCH --error=./LOG/merge3_%J.err
#SBATCH --output=./LOG/merge3_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.8.3_merge3.R
