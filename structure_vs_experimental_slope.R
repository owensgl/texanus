library(tidyverse)
library(stringr)
library(RColorBrewer)
ancestry_full <- data.frame(sample=as.character(),window=as.numeric(),window_start=as.numeric(),
                      window_end=as.numeric(),group=as.character(),contribution=as.numeric(),
                      percent=as.numeric(),
                      species=as.character(),tested=as.numeric(),chr=as.numeric())

for (n in 1:17){
  i <- str_pad(n, 2,pad="0")
  directory <- paste("/media/owens/Childs/RNAseq/texanus_linkage_structure/chr",i,sep="")
  file <- "results_ss.processed.ancestry"
  infofile <- "samplelist.txt"
  
  ancestry <- read.delim(paste(directory,file,sep="/"),header=T)
  info <- read.delim(paste(directory,infofile,sep="/"),header=F)
  colnames(info) <- c("name","species","tested")
  info$sample <- 1:nrow(info)
  ancestry <- merge(ancestry, info)
  ancestry %>% filter(tested == 0) -> ancestry
  ancestry$chr <- n
  ancestry_full <- rbind(ancestry_full, ancestry)
}

ancestry_summary <- ancestry_full %>% filter(species== 5) %>% 
  filter(group == 3) %>%
  group_by(chr, window_start,window_end,group,species) %>%
  #group_by(chr) %>%
  summarize(mean_percent=mean(percent)) %>%
  dplyr::rename(end = window_end,start=window_start) 

 abc_out_10seedbank_pointestimate_formerge <-  abc_out_10seedbank_pointestimate

 abc_out_10seedbank_pointestimate_formerge$chr <- as.integer(gsub("HanXRQChr", "", abc_out_10seedbank_pointestimate_formerge$chr))
 
pdf("texanus_RNAseq_structure_vs_relfitness.v0.pdf")
merge(ancestry_summary,
      abc_out_10seedbank_pointestimate_formerge %>% filter(parent == "deb")) %>% 
  ggplot(.,aes(x=mean_percent,y=rel_fitness_10seed, color=pop)) + geom_point() +
  geom_smooth(method=lm,se=F) + 
  scale_color_brewer(palette = "Set1",name="Location") +
  ylab("Relative_fitness") +
  xlab("Wild introgression percent") + 
  theme_bw()

brewer.pal(6, "Set1")
merge(ancestry_summary,
      abc_out_10seedbank_pointestimate_formerge %>% filter(parent == "deb")) %>%
  mutate(introgressed = if_else(mean_percent > 0.01, "yes","no")) %>% 
  ggplot(.,aes(x=pop, y=rel_fitness_10seed,fill=introgressed)) + geom_boxplot() + 
  theme_bw() + scale_fill_manual(values=brewer.pal(6, "Set1")[5:6])
dev.off()
summary(merge(ancestry_summary,
              abc_out_10seedbank_pointestimate_formerge %>% filter(parent == "deb")) %>%
  lm(mean_percent ~ rel_fitness_10seed + pop ,data=.))


pdf("texanus_RNAseq_structure_vs_relslope.v0.pdf")
merge(ancestry_summary
      ,windowrel.lm %>% filter(parent == "deb")) %>%  
  ggplot(.,aes(x=mean_percent,y=estimate, color=location)) + geom_point() +
  geom_smooth(method=lm,se=F) + 
  scale_color_brewer(palette = "Set1",name="Location") +
  ylab("Experimental selection slope") +
  xlab("Wild introgression percent") + 
  theme_bw()
dev.off()
summary(merge(ancestry_summary
              ,windowrel.lm %>% filter(parent == "deb")) %>%
          lm(estimate ~ mean_percent   ,data=.))


 
  
  
  
