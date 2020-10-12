#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --job-name=errF
#SBATCH --error=./LOG/errF_%J.err
#SBATCH --output=./LOG/errF_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.4_errF.R
