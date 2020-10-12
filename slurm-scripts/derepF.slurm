#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=350gb
#SBATCH --job-name=derepF
#SBATCH --error=./LOG/derepF_%J.err
#SBATCH --output=./LOG/derepF_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.2_derepF.R
