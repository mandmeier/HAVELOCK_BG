#================================================================================================
# SCRIPT 1 – filter_and_trim.R
#================================================================================================

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

### path to unzipped forward and reverse fastq files
miseq_path <- "/work/jyanglab/mmeier/HAVELOCK_BG/largedata/BG_16S"

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)
save(sampleNames, file = "cache/sampleNames.rda")

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)


### make 'filtered' folder and define the filenames for the filtered fastq.gz files.
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

save(filtFs, file = "cache/filtFs.rda")
save(filtRs, file = "cache/filtRs.rda")
print("LOG ==> saved filt")

#================================================================================================
### filter and trim fastq reads
#================================================================================================

### filter and trim using dada2
##### THIS TAKES ~1h per 100 f/r samples
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
  maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
  compress=TRUE, multithread=FALSE)

