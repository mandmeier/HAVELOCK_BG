#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=16gb
#SBATCH --job-name=mafttree
#SBATCH --error=./LOG/mafttree_%J.err
#SBATCH --output=./LOG/mafttree_%J.out

module load mafft/7.407
mafft cache/unaligned_ps_core5.fasta > cache/aligned_ps_core5.fasta

module load fasttree/2.1
FastTree -gtr -nt cache/aligned_ps_core5.fasta > cache/ps_core5.tre
