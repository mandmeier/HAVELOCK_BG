import sys
import os
import pandas as pd
import datetime
import time
import subprocess

## usage python /work/jyanglab/mmeier/HAVELOCK_BG/profiling/8.GWAS/run_gemma.py traits_150_stdN.txt


root_dir = "/work/jyanglab/mmeier/HAVELOCK_BG"

traits_file = sys.argv[1]
traits_path = os.path.join(root_dir,"cache/GWAS",traits_file)

snp_path = os.path.join(root_dir,"largedata/GWAS/282_genotype+Q+K/allchr_bisnp_n282_snpid_maf01_geno2")
kinship_path = os.path.join(root_dir,"cache/GWAS/Centered_IBS_kinship_GEMMA.txt")
pca_path = os.path.join(root_dir,"cache/GWAS/PCA.txt")


## make log directory
traits = os.path.splitext(traits_file)[0]
date_time = datetime.datetime.now().strftime("%y%m%d-%H%M%S")
log_path = os.path.join(root_dir,"LOG", f'{traits}_{date_time}')
os.mkdir(log_path)

outdir_path = os.path.join(root_dir,f'largedata/GWAS/out_{traits}_{date_time}')


def write_slurm(num):
    slurm_file = open("gemma.slurm", "w")
    slurm_file.write(
        f'#!/bin/sh\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=8\n#SBATCH --time=24:00:00\n#SBATCH --mem=16gb\n#SBATCH --job-name=T{num}\n#SBATCH --error={log_path}/T{num}_%J.err\n#SBATCH --output={log_path}/T{num}_%J.out\n\nmodule load python/3.8\n\ngemma-0.98 -bfile {snp_path} -k {kinship_path} -c {pca_path} -p {traits_path} -lmm 1 -n {num} -outdir {outdir_path} -o T{num} -miss 0.9 -r2 1 -hwe 0 -maf 0.01')
    slurm_file.close()


for trait in range(1,151):
    # check number of jobs
    jobs = subprocess.run(['squeue', '-u', 'mmeier'], stdout=subprocess.PIPE).stdout.decode('utf-8').count('\n')
    print(str(jobs) + "jobs")
    # if too many jobs wait 10 min
    while jobs > 980:
        time.sleep(600)
        jobs = subprocess.run(['squeue', '-u', 'mmeier'], stdout=subprocess.PIPE).stdout.decode('utf-8').count('\n')
        
    print("continue")   
    write_slurm(trait)
    os.system(f'echo "sbatch gemma.slurm for T{trait}"')
    os.system(f'sbatch gemma.slurm')


