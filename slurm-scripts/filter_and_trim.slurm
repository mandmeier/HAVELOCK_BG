#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=96gb
#SBATCH --job-name=filter
#SBATCH --error=./LOG/filter_%J.err
#SBATCH --output=./LOG/filter_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.1_filter_and_trim.R
