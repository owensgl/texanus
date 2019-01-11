#Plot the debilis alleles per sample

genotypes <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.15-35minor.parentagesites.dp10.tidygenos.txt",header=T,
                        colClasses = c("factor","factor","numeric","factor","numeric","factor","numeric"))

genotypes <- merge(genotypes, new_sites %>% select(sample,percent_new))
pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.15-35minor.parentagesites.dp10.tidygenos.pdf",height=4,width=6)
for (i in 1:4){
  genotypes %>% filter(location == poplist[i]) %>% select(gen) %>% unique() -> genlist
  for (j in 1:nrow(genlist)){
    for (k in 1:nrow(chrlist)){
      print(
      genotypes %>% filter(location == poplist[i],gen == genlist[j,]) %>%
        filter(!is.na(genotype)) %>%
        filter(chr == as.character(chrlist[k,])) %>%
        ggplot(.) + geom_point(aes(x=pos,y=as.factor(sample),col=as.factor(genotype))) +
        theme_minimal() +
        theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank()) +
        scale_colour_brewer(palette="Set1", 
                            name="Allelic state",
                            breaks=c("0", "1", "2"),
                            labels=c("0", "1", "2")) +
        ggtitle(paste(poplist[i], genlist[j,], chrlist[k,], sep=" ")) +
        facet_grid(~parent) 
      )
      }
  }
}
dev.off()




pdf("test.pdf",height=4,width=5)
genotypes %>% filter(location == "BFL",gen == "14") %>%
  filter(!is.na(genotype)) %>%
  filter(chr == "HanXRQChr02") %>%
  ggplot(.) + geom_point(aes(x=pos,y=as.factor(sample),col=as.factor(genotype))) +
  theme_minimal() +
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_colour_brewer(palette="Set1", 
                    name="Allelic state",
                    breaks=c("0", "1", "2"),
                    labels=c("0", "1", "2")) 
dev.off()
  
