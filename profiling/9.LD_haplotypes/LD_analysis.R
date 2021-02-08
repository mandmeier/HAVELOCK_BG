#LD plots

library("genetics")
library("LDheatmap")
library("tidyverse")


### find top bins

load("cache/top_taxa_bins_genes.rda")
unique(top_bins$tax_group)

hmp_bins <- top_bins %>%
  ungroup() %>%
  dplyr::select(chr, bin, sign_snps) %>%
  unique()

write_csv(hmp_bins, file = "cache/hmp_bins.csv")


data <- read.table("largedata/mdp_genotype.hmp.txt",head=T,com="",sep="\t",row=1,na.string="NN")

### find SNPs in 10kb window
ano <- filter(top_bins, tax_group == "Acinetobacter nosocomialis")

chr3 <- data %>%
  filter(chrom == 3) %>%
  filter(pos >= 186000000 & pos <= 189500000)


?LDheatmap

#gene <- data.frame(t(data[,11:ncol(data)]))

gene <- data.frame(t(chr3[,11:ncol(chr3)]))
gty <- makeGenotypes(gene,sep="")

rgb.palette <- colorRampPalette(rev(c("blue","orange" ,"red")), space = "rgb")
png(file="figures/LD_plots/LD_plot_chr3_1860-1895.png",res=300,width=1000,height=1000)
myld <- LDheatmap(gty,genetic.distances=chr3$pos,flip=TRUE, text=TRUE,color=rgb.palette(18))
dev.off()



myLDplot=function(file,LD=F,out=F){
  
  data=read.table(file,head=T,com="",sep="\t",row=1,na.string="NN")
  gene=data[,11:ncol(data)]
  gene=data.frame(t(gene))
  gty=makeGenotypes(gene,sep="")
  rgb.palette <- colorRampPalette(rev(c("blue","orange" ,"red")), space = "rgb")
  png(file=sub("txt","png",file),res=300,width=2000,height=2000)
  myld=LDheatmap(gty,genetic.distances=data$pos,flip=T,color=rgb.palette(18))
  dev.off()
  if(LD){
    write.table(myld$LDmatrix,file=sub(".txt","_LD.txt",file),sep="\t",quote=F)
  }
  if(out){ myld}else{NULL}
}
#example
myLDplot("my.hmp.txt",LD=T)