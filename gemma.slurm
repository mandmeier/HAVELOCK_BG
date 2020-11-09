#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=16gb
#SBATCH --job-name=T150
#SBATCH --error=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/traits_150_lowN_201108-094604/T150_%J.err
#SBATCH --output=/work/jyanglab/mmeier/HAVELOCK_BG/LOG/traits_150_lowN_201108-094604/T150_%J.out

module load python/3.8

gemma-0.98 -bfile /work/jyanglab/mmeier/HAVELOCK_BG/largedata/GWAS/282_genotype+Q+K/allchr_bisnp_n282_snpid_maf01_geno2 -k /work/jyanglab/mmeier/HAVELOCK_BG/cache/GWAS/Centered_IBS_kinship_GEMMA.txt -c /work/jyanglab/mmeier/HAVELOCK_BG/cache/GWAS/PCA.txt -p /work/jyanglab/mmeier/HAVELOCK_BG/cache/GWAS/traits_150_lowN.txt -lmm 1 -n 150 -outdir /work/jyanglab/mmeier/HAVELOCK_BG/largedata/GWAS/out_traits_150_lowN_201108-094604 -o T150 -miss 0.9 -r2 1 -hwe 0 -maf 0.01