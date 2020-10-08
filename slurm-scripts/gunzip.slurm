#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=48gb
#SBATCH --job-name=gunz
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/gunz_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/gunz_%J.out

gunzip /work/jyanglab/mmeier/HAVELOCK_BG/largedata/BG_16S/*.gz

