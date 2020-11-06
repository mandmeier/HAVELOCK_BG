# shorten.py

import sys
import pandas as pd

infile = sys.argv[1]
threshold = float(sys.argv[2])
outfile = sys.argv[3]

df = pd.read_csv(infile, sep="\t")
df = df[['rs','chr','ps','beta','se','p_wald']]
df = df[df.p_wald.le(threshold)]
df.to_csv(outfile, sep='\t', index=False)
