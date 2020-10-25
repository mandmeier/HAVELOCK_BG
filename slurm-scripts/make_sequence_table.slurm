#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=62gb
#SBATCH --job-name=table
#SBATCH --error=./LOG/table_%J.err
#SBATCH --output=./LOG/table_%J.out

module load R/3.5
Rscript profiling/1.16S_seq_processing/1.9_make_sequence_table.R


