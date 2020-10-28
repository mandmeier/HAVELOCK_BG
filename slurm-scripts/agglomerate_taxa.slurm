#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --job-name=glom
#SBATCH --error=./LOG/glom_%J.err
#SBATCH --output=./LOG/glom_%J.out

module load R/3.6
Rscript profiling/4.Find_core_taxa/agglomerate_taxa.R
