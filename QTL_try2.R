library(qtl)
library(tidyverse)
#This uses the scanone and scantwo functions to get qtl size and effect.
#deb <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.recode.csvr",estimate.map=FALSE)
deb <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.recode.imputed.csvr",estimate.map=FALSE)

deb.nt <- ntyped(deb, "mar")
todrop <- names(deb.nt[deb.nt < 20])
deb <- drop.markers(deb, todrop)

deb <- sim.geno(deb, step=0, n.draws=1000, err=0.01)
deb.ehk<- scanone(deb, pheno.col=1:25,method="ehk")
deb.perm <- scanone(deb, pheno.col=1:25,n.perm=1000,method="ehk")
deb.hk.2 <- scantwo(deb,pheno.col=1:25,method="hk")
deb.hk.2.perm <- scantwo(deb,pheno.col=1:25,method="hk",n.perm=1000)

plot(deb.ehk, lodcol=23:25,col=c("blue", "red","green"))
plot(deb.hk.2, chr=c(1,4,6,7,15))
summary(deb.perm, alpha=0.05)
allpeaks <- summary(deb.ehk, format="allpeaks", perms=deb.perm, alpha=0.05, pvalues=TRUE)
tidy_peaks <- data.frame(trait=as.character(),chr=as.character(),
                         left_edge=as.numeric(),right_edge=as.numeric(),
                         peak=as.numeric())
for (i in 1:25){
  pcol <- (i-1)*3+4
  chrcol <- 1
  loccol <- (i-1)*3+2
  trait <- colnames(allpeaks)[(i-1)*3+3]
  for (j in 1:nrow(allpeaks)){
    pvalue <- allpeaks[j, pcol]
    chr <- allpeaks[j,chrcol]
    loc <- allpeaks[j, loccol]
    if (pvalue < 0.05){
      interval <- bayesint(deb.ehk, chr, 0.95, lodcolumn=i)
      interval[3,2]
      tmp_peak <- data.frame(trait=as.character(trait),chr=as.character(chr),
                               left_edge=as.numeric(interval[1,2]),right_edge=as.numeric(interval[3,2]),
                               peak=as.numeric(loc))
      tidy_peaks <- rbind(tidy_peaks,tmp_peak)
      
    }
  }
}
tidy_peaks$chr <- factor(tidy_peaks$chr, levels= tidy_peaks$chr[order(as.numeric(as.character(tidy_peaks$chr)))])
levels(tidy_peaks$chr)
pdf("texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.bayesqtlpeaks.pdf",height=4,width=12)
tidy_peaks %>% 
  ggplot(.) + 
  geom_segment(aes(x=left_edge,xend=right_edge,y=as.factor(trait),
                                            yend=as.factor(trait),color=as.factor(trait))) +
  geom_point(aes(x=peak,y=as.factor(trait),color=as.factor(trait))) +
  facet_grid(~chr) +
  scale_y_discrete(name="", limits = rev(levels(tidy_peaks$trait))) +
  xlab("cM") + scale_color_discrete(name="Traits")
dev.off()
qtl <- makeqtl(deb, chr=c("05","06","13","15","17"), pos=c(48.9,19.9,10.5,44.7,56.8))
qtl <- makeqtl(deb, chr=c("05","06","13"), pos=c(48.9,19.9,10.5))

out.fq <- fitqtl(deb, qtl=qtl, formula=y~Q1+Q2+Q3)
summary(out.fq)
bayesint(deb.ehk, 13, 0.95, lodcolumn=2)



###Without imputation
library(qtl)
library(tidyverse)
#This uses the scanone and scantwo functions to get qtl size and effect.
debr <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.recode.csvr",estimate.map=FALSE)

debr.nt <- ntyped(debr, "mar")
todrop <- names(debr.nt[deb.nt < 20])
debr <- drop.markers(debr, todrop)

debr.ehk<- scanone(deb, pheno.col=1:25,method="ehk")
debr.perm <- scanone(deb, pheno.col=1:25,n.perm=1000,method="ehk")
summary(debr.ehk, format="allpeaks", perms=debr.perm, alpha=0.05, pvalues=TRUE)
qtl.r <- makeqtl(debr, chr=c("05","06","13","15","17"), pos=c(48.9,19.9,10.5,44.7,56.8))


