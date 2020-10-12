#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=72:00:00
#SBATCH --mem=500gb
#SBATCH --job-name=dadaR
#SBATCH --error=./LOG/dadaR_%J.err
#SBATCH --output=./LOG/dadaR_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.7_dadaR.R
