#For mapping contamination at the genome window level
library(tidyverse)
library(ggthemes)
library(gridExtra)

contam <- read_tsv("/home/owens/working/texanus/texanus.gatk.GBS2.BC1.LBJ.missing50.dp5.10mac.bi.contaminants.percent.txt",
                   col_types="cnncnnn")
#contam <- read_tsv("/home/owens/working/texanus/texanus.gatk.GBS2.BC1.LBJ.missing50.dp5.10mac.bi.contaminants.late10maf.all.percent.txt",
                   col_types="cnncnnn")


info <- read_tsv("/home/owens/working/texanus/texanus.sampleinfo.txt")

split_point <- 0.01
min_sites <- 40

hist.1 <- inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>%
  filter(percent_contam_homo < 0.25) %>%
  ggplot(.,aes(percent_contam_homo)) +geom_histogram(bins=50) +
  xlab("Proportion_homozygous_novel_sites") +
  geom_vline(xintercept = split_point,color="red")

hist.2 <- inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>%
  filter(percent_contam_het < 0.25) %>%
  ggplot(.,aes(percent_contam_het)) +geom_histogram(bins=50) +
  xlab("Proportion_heterozygous_novel_sites") +
  geom_vline(xintercept = split_point,color="red")

grid.arrange(hist.1,hist.2)


inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>% 
  filter(chr=="HanXRQChr01",window_start == 10000000) %>%
  ggplot(.) + geom_point(aes(x=percent_contam_het,y=percent_contam_homo))

  
inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>% 
  mutate(contaminated = case_when(percent_contam_homo > split_point ~ "2",
                                  percent_contam_het > split_point ~ "1",
                                  TRUE ~ "0")) %>% 
  mutate(ID = fct_reorder(ID,-gen)) %>% 
  ggplot(.) + geom_segment(aes(x=window_start/1000000,xend=window_end/1000000,y=ID,yend=ID,color=contaminated)) +
  facet_wrap(~chr,scales="free_x") +scale_color_brewer(palette = "Set1",name="Contamination") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("MB")

inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>%
  mutate(contaminated = case_when(percent_contam_homo > split_point ~ "2",
                                  percent_contam_het > split_point ~ "1",
                                  TRUE ~ "0")) %>% 
  group_by(ID, contaminated, gen) %>%
  tally() %>%
  spread(contaminated, n,fill=0) %>% 
  mutate(genome_contaminated = ( (`1` + (`2`*2)) / ((`0` + `1` + `2`)*2))) %>% 
  ggplot(.) + 
  geom_boxplot(aes(x=as.factor(gen),y=genome_contaminated*100),color="lightgrey", outlier.size=0) +
  geom_jitter(aes(x=as.factor(gen),y=genome_contaminated*100),width=0.2) +
  ylab("Percent genome from outside sources")


inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>%
  mutate(contaminated = case_when(percent_contam_homo > split_point ~ "2",
                                  percent_contam_het > split_point ~ "1",
                                  TRUE ~ "0")) %>% 
  group_by(ID, contaminated, gen) %>%
  tally() %>%
  spread(contaminated, n,fill=0) %>% 
  mutate(genome_contaminated = ( (`1` + (`2`*2)) / ((`0` + `1` + `2`)*2)),
         heterozygosity = (`1` / (`0` + `1` + `2`)) ) %>% 
  ggplot(.) + 
  facet_wrap(~gen) +
  coord_cartesian(ylim=c(0,1)) +
  geom_segment(aes(x=0,xend=0.5,y=0,yend=1),color="grey") +
  geom_segment(aes(x=0.5,xend=1,y=1,yend=0),color="grey") +
  geom_segment(aes(x=0,xend=1,y=0,yend=0),color="grey") +
  theme_few() +
  geom_point(aes(x=genome_contaminated,y=heterozygosity)) +
  xlab("Contamination percent") + ylab("Contamination Heterozygosity")



inner_join(contam,info) %>% 
  filter(total_sites > min_sites) %>% 
  mutate(contaminated = case_when(percent_contam_homo > split_point ~ 2,
                                  percent_contam_het > split_point ~ 1,
                                  TRUE ~ 0)) %>% 
  group_by(chr,gen, window_start, window_end) %>%
  summarize(mean_contam = mean(contaminated,na.rm=T),count=n()) %>% 
  filter(count > 10) %>%
  mutate(window_middle = window_start+(window_end - window_start)/2) %>% 
  ggplot(.,aes(x=window_middle/1000000,y=mean_contam,group=gen,color=gen)) + geom_line() +
  facet_wrap(~chr,scales="free_x",nrow=3) + scale_color_distiller(palette = "Spectral",name="Generation") +
  theme_few() + xlab("MB") +ylab("Mean contamination")
  





