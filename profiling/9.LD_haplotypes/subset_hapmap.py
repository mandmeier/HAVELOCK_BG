import pandas as pd

hmp_bins = 'cache/hmp_bins.csv'
#hmp_bins = 'temp/hmp_bins_test.csv'
hmp_large = 'largedata/hmp321_282_agpv4_maf005_miss03.hmp.txt'
#hmp_large = 'largedata/hmp.txt'
hmp_out = 'cache/top_snps_hmp.txt'

# get bins from file
print('get bins from file...')
bins = pd.read_csv(hmp_bins)
bins['index'] = 'chr' + bins.chr.astype(str) + '_' + (bins.bin - 1 ).astype(str)

top_snps_hmp = pd.DataFrame()

# read data
for df in pd.read_csv(hmp_large,sep='\t', chunksize=100000):
    df['index'] = 'chr' + df.chrom.astype(str) + '_' + df.pos.astype(str).map(lambda x: x[:-4])
    # get top snps from bins
    top_snps = df[df['index'].isin(bins['index'])]
    print(f'found {top_snps.shape[0]} SNPs in {df.shape[0]} rows')
    #print(top_snps)
    top_snps_hmp = top_snps_hmp.append(top_snps)

# save to tab delimited file
print('save txt file...')
top_snps_hmp.to_csv(hmp_out, sep='\t', index=False)