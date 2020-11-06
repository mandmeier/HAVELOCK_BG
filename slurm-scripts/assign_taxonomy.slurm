#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --job-name=tax
#SBATCH --error=./LOG/tax_%J.err
#SBATCH --output=./LOG/tax_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.10_assign_taxonomy.R
