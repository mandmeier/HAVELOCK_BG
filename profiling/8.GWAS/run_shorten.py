### shorten huge GWAS output files, retain only highly significant snps
## usage
## python run_shorten.py out_traits_150_lowN_201104-185415 5

import sys
import os
import subprocess

root_dir = "/work/jyanglab/mmeier/HAVELOCK_BG"
largedat_dir = sys.argv[1]
thresharg = sys.argv[2]
threshold = float(f'1E-{thresharg}')

#largedat_dir = "out_traits_stdN_200813-192612"

## input direcory
largedat_path = os.path.join(root_dir,"largedata/GWAS",largedat_dir)

## make log directory

log_path = os.path.join(root_dir,"LOG", f'shorten_{thresharg}_{largedat_dir}')

if not os.path.isdir(log_path):
    os.mkdir(log_path)


## make output directory
shortdat_path = os.path.join(root_dir, "largedata/GWAS", f'short_{thresharg}_{largedat_dir}')
if not os.path.isdir(shortdat_path):
    os.mkdir(shortdat_path)



def write_slurm(trait, item, threshold, outfile):
    slurm_file = open("shorten.slurm", "w")
    slurm_file.write(
        f'#!/bin/sh\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=8\n#SBATCH --time=12:00:00\n#SBATCH --mem=16gb\n#SBATCH --job-name={trait}\n#SBATCH --error={log_path}/{trait}_%J.err\n#SBATCH --output={log_path}/{trait}_%J.out\n\nmodule load python/3.8\n\npython /work/jyanglab/mmeier/HAVELOCK_BG/profiling/8.GWAS/shorten.py {infile} {threshold} {outfile}')
    slurm_file.close()


for item in os.listdir(largedat_path):
    if item.endswith('.assoc.txt'):

        infile = os.path.join(largedat_path, item)

        # check number of jobs
        jobs = subprocess.run(['squeue', '-u', 'mmeier'], stdout=subprocess.PIPE).stdout.decode('utf-8').count('\n')
        print(str(jobs) + "jobs")
        # if too many jobs wait 10 min
        while jobs > 990:
            time.sleep(600)
            jobs = subprocess.run(['squeue', '-u', 'mmeier'], stdout=subprocess.PIPE).stdout.decode('utf-8').count('\n')

        trait = item.split('.assoc')[0]

        outfile = os.path.join(shortdat_path,f'{trait}.assoc.txt')

        write_slurm(trait, infile, threshold, outfile)
        os.system(f'echo "sbatch shorten.slurm for {trait}"')
        os.system(f'sbatch -p jyanglab --licenses=common --ntasks=1 --mem 16G --time=12:00:00 shorten.slurm')
        print(trait)


