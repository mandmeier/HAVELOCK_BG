#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=8:00:00
#SBATCH --mem=64gb
#SBATCH --job-name=m64
#SBATCH --error=./LOG/m64_%J.err
#SBATCH --output=./LOG/m64_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.8.1_merge1.R
