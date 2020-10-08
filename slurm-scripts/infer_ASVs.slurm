#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=8gb
#SBATCH --job-name=ASVs
#SBATCH --error=./LOG/ASVs_%J.err
#SBATCH --output=./LOG/ASVs_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.2_infer_ASVs.R

