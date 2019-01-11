library(qtl)
library(qtlcharts)
iplotScanone(deb.em, deb,fillgenoArgs=list(method="argmax", error.prob=0.01))
# deb <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.recode.csvr",estimate.map=FALSE)
# ann1 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann1sites.recode.csvr",estimate.map=FALSE)
# ann2 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann2sites.recode.csvr",estimate.map=FALSE)
# ann3 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann3sites.recode.csvr",estimate.map=FALSE)
deb <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2debsites.recode.combined.csvr",estimate.map=FALSE)
ann2 <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.gatk.GBS1.BC1.maf10.bi.GBS2ann2sites.recode.combined.csvr",estimate.map=FALSE)
snpchip <- read.cross(format="csvr", "/home/owens/working/texanus/","texanus.snplist.data.cm.csvr",estimate.map=FALSE)

deb.nt <- ntyped(deb, "mar")
todrop <- names(deb.nt[deb.nt < 20])
deb <- drop.markers(deb, todrop)


deb <- calc.errorlod(deb, error.prob=0.01)
print(toperr <- top.errorlod(deb, cutoff=6))

deb.clean <- deb
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  deb.clean$geno[[chr]]$data[deb$pheno$id==id, mar] <- NA
}


ann2.nt <- ntyped(ann2, "mar")
todrop <- names(ann2.nt[ann2.nt < 20])
ann2 <- drop.markers(ann2, todrop)


ann2 <- calc.errorlod(ann2, error.prob=0.01)
print(toperr <- top.errorlod(ann2, cutoff=6))

ann2.clean <- ann2
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  ann2.clean$geno[[chr]]$data[ann2$pheno$id==id, mar] <- NA
}


deb <- calc.genoprob(deb,step=0)
ann2 <- calc.genoprob(ann2,step=0)


nperm=500
pdf("texanus.gatk.GBS1.BC1.maf10.bi.GBS2.QTLs.combined.v5.pdf", height=10,width=10)
par(mfrow=c(3,1))
for (i in 1:25){
  deb.em <- scanone(deb, pheno.col=i)
  deb.perm <- scanone(deb, pheno.col=i,n.perm=nperm)
  plot(deb.em)
  title(paste("Debilis -",colnames(tmp1$pheno[i])))
  abline(h=summary(deb.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(deb.perm)[[2]],col="orange",lty=2,lwd=2)
  
  snpchip.em <- scanone(snpchip, pheno.col=i)
  snpchip.perm <- scanone(snpchip, pheno.col=i,n.perm=nperm)
  plot(snpchip.em)
  abline(h=summary(snpchip.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(snpchip.perm)[[2]],col="orange",lty=2,lwd=2)
  title(paste("SNPchip -",colnames(tmp1$pheno[i])))
  
  # ann1.em <- scanone(ann1, pheno.col=i)
  # ann1.perm <- scanone(ann1, pheno.col=i,n.perm=nperm)
  # plot(ann1.em)
  # title(paste("Ann1 -",colnames(tmp1$pheno[i])))
  # abline(h=summary(ann1.perm)[[1]],col="green",lty=2,lwd=2)
  # abline(h=summary(ann1.perm)[[2]],col="orange",lty=2,lwd=2)
  
  ann2.em <- scanone(ann2, pheno.col=i)
  ann2.perm <- scanone(ann2, pheno.col=i,n.perm=nperm)
  
  plot(ann2.em)
  title(paste("Ann2 -",colnames(tmp1$pheno[i])))
  abline(h=summary(ann2.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(ann2.perm)[[2]],col="orange",lty=2,lwd=2)
  
  # ann3.em <- scanone(ann3, pheno.col=i)
  # ann3.perm <- scanone(ann3, pheno.col=i,n.perm=nperm)
  # plot(ann3.em)
  # abline(h=summary(ann3.perm)[[1]],col="green",lty=2,lwd=2)
  # abline(h=summary(ann3.perm)[[2]],col="orange",lty=2,lwd=2)
  # title(paste("Ann3 -",colnames(tmp1$pheno[i])))
  # plot(1, 1, type = "n")

  
}
dev.off()



deb.clean <- sim.geno(deb.clean, step=0, n.draws=1000, err=0.01)
ann2.clean <- sim.geno(ann2.clean, step=0, n.draws=1000, err=0.01)

nperm=1000

pdf("texanus.gatk.GBS1.BC1.maf10.bi.GBS2.QTLsEffect.combined.v5.pdf", height=10,width=10)
par(mfrow=c(2,2))
for (i in 1:25){
  deb.em <- scanone(deb.clean, pheno.col=i)
  deb.perm <- scanone(deb.clean, pheno.col=i,n.perm=nperm)
  
  plot(deb.em)
  title(paste("Debilis -",colnames(tmp1$pheno[i])))
  abline(h=summary(deb.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(deb.perm)[[2]],col="orange",lty=2,lwd=2)
  
  effectscan(deb.clean,pheno.col=i,get.se=T)
  title(paste("Debilis -",colnames(tmp1$pheno[i]),"Effect Size"))
  
  ann2.em <- scanone(ann2.clean, pheno.col=i)
  ann2.perm <- scanone(ann2.clean, pheno.col=i,n.perm=nperm)
  
  plot(ann2.em)
  title(paste("Annuus BC parent -",colnames(tmp1$pheno[i])))
  abline(h=summary(ann2.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(ann2.perm)[[2]],col="orange",lty=2,lwd=2)
  
  effectscan(ann2.clean,pheno.col=i,get.se=T)
  title(paste("Annuus BC parent -",colnames(tmp1$pheno[i]),"Effect Size"))
} 
dev.off()
  
scanone(ann1,pheno.col = 1:25)

