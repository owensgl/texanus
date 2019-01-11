library(qtl)
library(qtlcharts)
iplotScanone(deb.em, deb,fillgenoArgs=list(method="argmax", error.prob=0.01))
deb <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.recode.csvr",estimate.map=FALSE)
ann1 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann1sites.recode.csvr",estimate.map=FALSE)
ann2 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann2sites.recode.csvr",estimate.map=FALSE)
ann3 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann3sites.recode.csvr",estimate.map=FALSE)
snpchip <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.snplist.data.cm.csvr",estimate.map=FALSE)

deb.nt <- ntyped(deb, "mar")
todrop <- names(deb.nt[deb.nt < 20])
deb <- drop.markers(deb, todrop)

ann1.nt <- ntyped(ann1, "mar")
todrop <- names(ann1.nt[ann1.nt < 20])
ann1 <- drop.markers(ann1, todrop)

ann2.nt <- ntyped(ann2, "mar")
todrop <- names(ann2.nt[ann2.nt < 20])
ann2 <- drop.markers(ann2, todrop)
ann3.nt <- ntyped(ann3, "mar")
todrop <- names(ann3.nt[ann3.nt < 20])
ann3 <- drop.markers(ann3, todrop)

deb <- calc.genoprob(deb,step=0)
ann1 <- calc.genoprob(ann1,step=0)
ann2 <- calc.genoprob(ann2,step=0)
ann3 <- calc.genoprob(ann3,step=0)
snpchip <- calc.genoprob(snpchip,step=0)

snpchip.em <- scanone(snpchip, pheno.col=i)
plot(snpchip.em)
plot(deb.em)

pdf("texanus.gatk.GBS1.BC1.maf10.bi.GBS2.QTLs.v4.pdf", height=10,width=10)
par(mfrow=c(3,2))
for (i in 1:25){
  deb.em <- scanone(deb, pheno.col=i)
  deb.perm <- scanone(deb, pheno.col=i,n.perm=1000)
  plot(deb.em)
  title(paste("Debilis -",colnames(tmp1$pheno[i])))
  abline(h=summary(deb.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(deb.perm)[[2]],col="orange",lty=2,lwd=2)
  
  snpchip.em <- scanone(snpchip, pheno.col=i)
  snpchip.perm <- scanone(snpchip, pheno.col=i,n.perm=1000)
  plot(snpchip.em)
  abline(h=summary(snpchip.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(snpchip.perm)[[2]],col="orange",lty=2,lwd=2)
  title(paste("SNPchip -",colnames(tmp1$pheno[i])))
  
  ann1.em <- scanone(ann1, pheno.col=i)
  ann1.perm <- scanone(ann1, pheno.col=i,n.perm=1000)
  plot(ann1.em)
  title(paste("Ann1 -",colnames(tmp1$pheno[i])))
  abline(h=summary(ann1.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(ann1.perm)[[2]],col="orange",lty=2,lwd=2)
  
  ann2.em <- scanone(ann2, pheno.col=i)
  ann2.perm <- scanone(ann2, pheno.col=i,n.perm=1000)
  
  plot(ann2.em)
  title(paste("Ann2 -",colnames(tmp1$pheno[i])))
  abline(h=summary(ann2.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(ann2.perm)[[2]],col="orange",lty=2,lwd=2)
  
  ann3.em <- scanone(ann3, pheno.col=i)
  ann3.perm <- scanone(ann3, pheno.col=i,n.perm=1000)
  plot(ann3.em)
  abline(h=summary(ann3.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(ann3.perm)[[2]],col="orange",lty=2,lwd=2)
  title(paste("Ann3 -",colnames(tmp1$pheno[i])))
  plot(1, 1, type = "n")

  
}
dev.off()

max(snpchip.em)
mar <- find.marker(snpchip, chr="13", pos=51.9)
plotPXG(snpchip, marker=mar,pheno.col=i,infer=F)

deb.tmp <- scanone(deb, pheno.col=1,n.perm=1000)
plot(deb.tmp)
summary(deb.tmp)  #print the thresholds
plot(deb.em)
abline(h=summary(deb.tmp)[[1]],col="green",lty=2,lwd=2)
abline(h=summary(deb.tmp)[[2]],col="orange",lty=2,lwd=2)

pdf("texanus.gatk.GBS1.BC1.maf10.bi.GBS2.QTLsEffect.v4.pdf", height=10,width=10)
par(mfrow=c(2,1))
for (i in 1:25){
  deb.em <- scanone(deb, pheno.col=i)
  
  plot(deb.em)
  title(paste("Debilis -",colnames(tmp1$pheno[i])))
  
  effectscan(deb,pheno.col=i,get.se=T)
  title(paste("Debilis -",colnames(tmp1$pheno[i]),"Effect Size"))
} 
dev.off()
  
scanone(ann1,pheno.col = 1:25)

