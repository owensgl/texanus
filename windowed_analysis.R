
#window <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.10mb.txt",header=T,
#                        colClasses = c("factor","factor","factor","numeric","numeric","factor","numeric")))
window <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.txt",header=T,
                     colClasses = c("character","character","character","numeric","numeric","character","character"))


info <- read.delim("/home/owens/working/texanus/texanus.sampleinfo.txt",header=T)
info %>% select(ID, gen, location) -> info
colnames(info) <- c("sample","gen","location")
window <- merge(window, info)
bw_pallet <- c("#cccccc","#808080","#0c0c0c","#ffffff")

pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.bw.pdf",height=4,width=6)
for (i in 1:4){
  window %>% filter(location == poplist[i]) %>% select(gen) %>% unique() -> genlist
  for (j in 1:nrow(genlist)){
    for (k in 1:nrow(chrlist)){
      print(
        window %>%
          filter(location == poplist[i],gen == genlist[j,]) %>% 
          filter(chr == as.character(chrlist[k,])) %>%
          ggplot(.) +
          geom_segment(aes(x=start, xend=end,y=as.factor(sample),yend=as.factor(sample),color=call),size=2) +
          facet_grid(~parent) +
          scale_colour_manual(values=bw_pallet, 
                              name="Allelic state",
                              breaks=c("0", "1", "2","Draw"),
                              labels=c("0", "1", "2","Unknown")) +
          theme_minimal() +
          theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()) +
          ggtitle(paste(poplist[i], genlist[j,], chrlist[k,],"10MB window", sep=" "))
          
      )
    }
  }
}
dev.off()

pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.haplotypesobserved.pdf")
for (i in 1:4){
  print(
    window %>% filter(call != "Draw") %>%
      filter(location == poplist[i]) %>%
      group_by(sample,chr,start, end) %>%  mutate(count =n()) %>%
      filter(count == 4) %>% droplevels() %>%
      mutate(sum = sum(as.numeric(as.character(call)))) %>%
      ungroup() %>%
      group_by(sample) %>%
      mutate( percent_2 = sum(sum[sum == 2])/n()) %>% 
      ungroup() %>%
      arrange( desc(percent_2)) %>% 
      mutate(sample = factor(sample, levels = as.character(sample))) %>%
      ggplot(.) +
      geom_histogram(aes(x=sample,fill=as.factor(sum)),position="fill",stat="count") + 
      facet_wrap(~gen,drop=T,scales = "free_x") + 
      scale_fill_brewer(palette="Set1", 
                        name="Observed haplotypes") +
      theme(axis.text.x=element_blank()) +
      ylab("Proportion of genome with X observed haplotypes") +
      xlab("Generation") + 
      ggtitle(poplist[i])
  )
}
dev.off()  






pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.haplotypesonchr.pdf",height=4,width=6)
for (i in 1:4){
  genotypes %>% filter(location == poplist[i]) %>% select(gen) %>% unique() -> genlist
  for (j in 1:nrow(genlist)){
    for (k in 1:nrow(chrlist)){
      print(
        window %>% filter(call != "Draw") %>%
          filter(location == poplist[i],gen == genlist[j,]) %>% 
          filter(chr == as.character(chrlist[k,])) %>%
          group_by(sample,chr,start, end) %>%  mutate(count =n()) %>%
          filter(count == 4) %>% mutate(sum = sum(as.numeric(as.character(call))),
                                        callnumeric = as.numeric(as.character(call))) %>%
          ggplot(.) + 
          geom_segment(aes(x=start, xend=end,y=as.factor(sample),yend=as.factor(sample),color=as.factor(sum)),size=2) +
          scale_color_brewer(palette="Set1", 
                             name="Observed haplotypes") +
          theme_minimal() +
          theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()) +
          ggtitle(paste(poplist[i], genlist[j,], chrlist[k,],"10MB window haplotypecount", sep=" "))
      )
    }
  }
}
dev.off()

window.summarized <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.10mb.summarized.txt",header=T)

window.summarized <- merge(window.summarized,info)
new_sites <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.newalleles.txt",header=T,
                        colClasses = c("factor","factor","numeric","numeric","numeric","numeric"))
window.summarized <- merge(window.summarized,new_sites) 

pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.10mb.contamination.pdf")
window.summarized %>% ggplot(.) + 
  geom_jitter(aes(x=gen,y=percent_missing)) + facet_wrap(~location) +
  ylab("Percent of haplotypes missing") + xlab("Generation") +
  ggtitle("Percent contamination based on haplotype counts") + 
  theme_minimal()

window.summarized %>% ggplot(.) + 
  geom_jitter(aes(x=percent_new,y=percent_missing,col=location)) +
  scale_color_brewer(palette="Set1", name="Location") +
  theme_minimal() + 
  xlab("Proportion novel alleles") + 
  ylab("Percent of haplotypes missing") +
  ggtitle("Comparing methods to detect contamination")

dev.off()
  


