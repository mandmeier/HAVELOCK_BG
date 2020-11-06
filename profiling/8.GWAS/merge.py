# merge shortened GWAS output for all traits into one file

# usage:
# python merge.py <directory with individual GWAS output files>

import sys
import os, glob
import pandas as pd

arg = sys.argv[1]

root_dir = "/work/jyanglab/mmeier/HAVELOCK_BG"

assoc_dir = os.path.join(root_dir, "largedata/GWAS",  arg)

all_files = glob.glob(os.path.join(assoc_dir, "*assoc.txt"))

fname = "merged_" + arg + ".csv"

outfile = os.path.join(root_dir, "largedata/GWAS/", fname)

all_df = []
for f in all_files:
    print(f)
    df = pd.read_csv(f, sep='\t')
    df['trait'] = f.split('/')[-1].split('.')[0]
    all_df.append(df)

merged_df = pd.concat(all_df, ignore_index=True, sort=True)
merged_df.to_csv(outfile)

print(outfile)
