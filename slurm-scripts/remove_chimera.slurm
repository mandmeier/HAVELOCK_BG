#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=62gb
#SBATCH --job-name=chimera
#SBATCH --error=./LOG/chimera_%J.err
#SBATCH --output=./LOG/chimera_%J.out

module load R/3.6
Rscript profiling/1.16S_seq_processing/1.10_remove_chimera.R




