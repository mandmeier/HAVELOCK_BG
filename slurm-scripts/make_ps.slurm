#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --job-name=make_ps
#SBATCH --error=./LOG/make_ps_%J.err
#SBATCH --output=./LOG/make_ps_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.11_make_ps.R
