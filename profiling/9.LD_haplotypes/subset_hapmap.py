import pandas as pd

hmp_large = '/common/jyanglab/shared/Gen_Xu/282_genotype+Q+K/hmp321_282_agpv4_maf005_miss03.hmp.txt'
hmp_bins = 'cache/hmp_bins.csv'
hmp_out = 'cache/top_snps_hmp.txt'

# read data
df = pd.read_csv(hmp_large, sep="\t")
df['index'] = 'chr' + df.chrom.astype(str) + '_' + df.pos.astype(str).map(lambda x: x[:-4])

# get bins from file
bins = pd.read_csv(hmp_bins)
bins['index'] = 'chr' + bins.chr.astype(str) + '_' + (bins.bin - 1 ).astype(str)

top_snps_hmp = df[df['index'].isin(bins['index'])]

#save to tab delimited file
top_snps_hmp.to_csv(hmp_out, sep='\t', index=False)

