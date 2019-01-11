library(qtl)
chrlist <- read.delim("/home/owens/working/texanus/rqtl/chrlist.txt",header=F)

wgs_snps <- read.delim("/home/owens/working/texanus/WGS/texanus.gatk.GBS1.BC1moredata.maf15-35.bi.dp10.missing25.WGScompared.txt",header=T,comment.char="#")
wgs_snps %>% select(chr, pos) -> wgs_snps_pos
wgs_snps_pos$type <- "Deb"
marker_identity <- data.frame(chr=as.character(),pos=as.numeric(),parent=as.character())
full_marker_set <- data.frame(chr=as.character(),pos=as.numeric(),lg=as.numeric())
pdf("texanus.gatk.GBS1.BC1moredata.maf15-35.bi.dp10.missing25.linkageplots.pdf")
for (i in 1:nrow(chrlist)){
  mapthis <- read.cross(format="csvr", "/home/owens/working/texanus/rqtl",paste("texanus.gatk.GBS1.BC1moredata.maf15-35.bi.dp10.missing25.",chrlist[i,],".csvr",sep=""),estimate.map=FALSE,sep="\t")
  nmarkers <- summary.map(mapthis)$n.mar[1]
  min.markers <- round(nmarkers * .75)
  mapthis <- subset(mapthis, ind=(ntyped(mapthis)>min.markers))
  nt.bymar <- ntyped(mapthis, "mar")
  todrop <- names(nt.bymar[nt.bymar < 400])
  mapthis <- drop.markers(mapthis, todrop)
  mapthis <- est.rf(mapthis)
  mapthis <- formLinkageGroups(mapthis, max.rf=0.10, min.lod=15, reorgMarkers=TRUE)
  plotRF(mapthis, chr=c(1:10),alternate.chrid=TRUE, main=paste(chrlist[i,]))
  lg <- data.frame(chr=as.character(),pos=as.numeric(),lg=as.numeric())
  for (j in 1:10){
    tmp <- markernames(mapthis, chr=j)
    tmp <- as.data.frame(tmp)
    tmp %>% separate(tmp, c("chr","pos")) -> tmp
    tmp$lg <- j
    lg <- rbind(lg,tmp)
  }
  lg <- merge(lg, wgs_snps_pos, all.x=TRUE)
  lg[which(is.na(lg$type)),]$type <- "Ann"
  min_deb <- 0.1
  print(
    lg %>% ggplot(.) + geom_bar(aes(x=as.factor(lg),fill=type),position = "fill") +
      ylab("Proportion Debilis sites") +
      xlab("Linkage group") +
      ggtitle(paste(chrlist[i,])) +
      geom_hline(yintercept=min_deb)
  )
  all_markers <- lg %>% select(chr, pos, lg)
  full_marker_set <- rbind(full_marker_set,all_markers)
}
dev.off()