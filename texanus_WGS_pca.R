library(SNPRelate)
library(tidyverse)
library(stringr)

vcf_filename<- c("/home/owens/working/texanus/WGS/texanusWGS.gatk.vcf.gz")
gds_filename<- c("/home/owens/working/texanus/WGS/texanusWGS.gatk.gds")
#Convert your vcf to gds for use with snprelate
snpgdsVCF2GDS(vcf_filename, gds_filename,  method="biallelic.only",ignore.chr.prefix="HanXRQChr")

genofile <- snpgdsOpen(gds_filename)
#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

snpset.id <- unlist(snpset_pruned)
#Run the PCA
pca <- snpgdsPCA(genofile, num.thread = 2, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.15, maf = 0.05,autosome.only = F)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

sampleinfo <- read.delim("/home/owens/working/texanus/WGS/WGS.poplist.txt",header=F)
colnames(sampleinfo) <- c("name","species")
#Make a dataframe of your PCA results
tab <- data.frame(name = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

#Merge the sampleinfo into that
tab <- merge(tab, sampleinfo)

#Plot a PCA image
ggplot(data=tab,aes(EV1,EV2)) + geom_point(aes(color=species))




####Run chr by chr
prefix<- "/home/owens/working/texanus/WGS/texanusWGS.gatk.chr"


for (i in 1:17){
  n = str_pad(i, 2, pad = "0")
  snpgdsVCF2GDS(paste(prefix, n, ".vcf.gz",sep=""), paste(prefix,n,".gds",sep=""),  method="biallelic.only",ignore.chr.prefix="HanXRQChr")
  genofile <- snpgdsOpen(paste(prefix,n,".gds",sep=""),allow.duplicate=T)
  #Prune for linkage
  snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)
  
  snpset.id <- unlist(snpset_pruned)
  #Run the PCA
  pca <- snpgdsPCA(genofile, num.thread = 2, eigen.cnt = 16, missing.rate = 0.15, maf = 0.05,autosome.only = F)
  pc.percent <- pca$varprop*100
  head(round(pc.percent, 2))
  
  sampleinfo <- read.delim("/home/owens/working/texanus/WGS/WGS.poplist.txt",header=F)
  colnames(sampleinfo) <- c("name","species")
  #Make a dataframe of your PCA results
  tab <- data.frame(name = pca$sample.id,
                    EV1 = pca$eigenvect[,1],    # the first eigenvector
                    EV2 = pca$eigenvect[,2],    # the second eigenvector
                    EV3 = pca$eigenvect[,3],
                    EV4 = pca$eigenvect[,4],
                    stringsAsFactors = FALSE)
  
  #Merge the sampleinfo into that
  tab <- merge(tab, sampleinfo)
  
  #Plot a PCA image
  ggplot(data=tab,aes(EV1,EV2)) + geom_point(aes(color=species))
}



