#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=32gb
#SBATCH --job-name=H2
#SBATCH --error=./LOG/H2_%J.err
#SBATCH --output=./LOG/H2_%J.out

module load R/3.5
Rscript profiling/6.Heritability_analysis/asv_level_H2.R
