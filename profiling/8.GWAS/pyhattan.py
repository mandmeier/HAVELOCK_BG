import os
import sys
from Pyhattan import FormatData, GenerateManhattan

arg = sys.argv[1]
cwd = os.getcwd()
input_dir = os.path.join(cwd, arg)

output_dir = os.path.join(cwd, 'plots')
os.mkdir(output_dir)

for item in os.listdir(input_dir):
    if item.endswith('.assoc.short.txt'):
        trait = item.split('.assoc')[0]
        print(trait)
        infile = os.path.join(input_dir, item)
        #print(infile)
        
        outfile = os.path.join(output_dir,f'{trait}_man.png')
        #print(outfile)

        dat_gemma_output = FormatData(infile)
        GenerateManhattan(dat_gemma_output, export_path=outfile,significance=5, colors=['#00008B', '#006400'])


