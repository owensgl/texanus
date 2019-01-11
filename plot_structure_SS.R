library(tidyverse)
library(stringr)

results <- data.frame(sample=as.character(),loci=as.numeric(),bp_start=as.numeric(),
                      bp_end=as.numeric(),cm_start=as.numeric(),cm_end=as.numeric(),
                      group=as.character(),y_start=as.numeric(),y_end=as.numeric(),
                      species=as.character(),tested=as.numeric(),chr=as.numeric())
for (n in 1:17){
  i <- str_pad(n, 2,pad="0")
  directory <- paste("/media/owens/Childs1/RNAseq/texanus_linkage_structure/chr",i,sep="")
  file <- "results_ss.processed.rect"
  infofile <- "samplelist.txt"
  
  stru <- read.delim(paste(directory,file,sep="/"),header=T)
  info <- read.delim(paste(directory,infofile,sep="/"),header=F)
  colnames(info) <- c("name","species","tested")
  info$sample <- 1:nrow(info)
  stru <- merge(stru, info)
  stru %>% filter(tested == 0) -> stru
  stru$chr <- n
  results <- rbind(results, stru)
}

pdf("texanus_RNAseq_structure.v0.pdf")
for (n in 1:17){
  print(
  results %>% filter(tested == 0,species == 5) %>%
    filter(!is.na(cm_start)) %>%
    filter(chr == n) %>%
    ggplot(.) + geom_rect(aes(xmin=cm_start,xmax=cm_end,ymin=y_start,ymax=y_end,
                              fill=as.factor(group))) +
    scale_fill_brewer(palette="Set1",name="Species",
                      labels=c("H. annuus", "H. argophyllus", "H. debilis", "H. petiolaris")) + 
    facet_grid(name~.) +
    ggtitle(paste("Chromosome",n,sep=" "))
  )
}
dev.off()

  
