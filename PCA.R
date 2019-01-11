library(SNPRelate)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(igraph)
folder <- "/home/owens/working/texanus"
gds.file <- "texanus.gatk.noQTL.800misscount.5dp.10mac.gds"
vcf.file <- "texanus.gatk.noQTL.800misscount.5dp.10mac.vcf"

snpgdsVCF2GDS(paste(folder,vcf.file,sep="/"), paste(folder,gds.file,sep="/"), method="biallelic.only", ignore.chr.prefix = "HanXRQChr")

genofile <- snpgdsOpen(paste(folder,gds.file,sep="/"))
snpgdsSummary(paste(folder,gds.file,sep="/"))
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=5)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


sampleid <- read.delim(paste(folder,"/texanus.noQTL.samplelist.txt",sep=""),header=F)
colnames(sampleid) <- "ID"
popinfo <- read.delim(paste(folder,"/texanus.sampleinfo.txt",sep=""),header=T)

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

tmp.tab <- tab[,c(1,2,3)]
colnames(tmp.tab) <- c("Taxa","PC1","PC2")
tab <- cbind(tab, sampleid)
tab2 <- merge(tab,popinfo)



ggplot(tab2 %>% filter(gen >= 1)) + geom_point(aes(x=EV1,y=EV2,color=gen),alpha=0.5) + facet_wrap(~location) +
  theme_bw() + scale_colour_distiller(palette="Spectral") 

  
  
snpgdsClose(genofile)
