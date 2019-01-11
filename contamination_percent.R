library(ggthemes)
library(tidyverse)
library(gridExtra)
new_sites <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.newalleles.txt",header=T,
                        colClasses = c("factor","factor","numeric","numeric","numeric","numeric"))
pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.newalleles.pdf")
new_sites %>% filter(gen != "1") %>%
  ggplot(.) + geom_jitter(aes(x=gen,y=percent_new)) +
  facet_wrap(~location) + theme_few() +
  ylab("Percent of variants not found in generation 1") + 
  xlab("Generation") +
  ggtitle("Contamination from outside sources")
new_sites %>% filter(gen != "1") %>% 
  ggplot(.) + geom_bar(aes(as.factor(gen),fill=contaminated),stat_count=T) +
  facet_wrap(~location) + theme_few() + scale_fill_brewer(palette="Set1") +
  xlab("Generation")+ ylab("Samples")
dev.off()

new_sites %>% mutate(contaminated = ifelse(percent_new > 0.005, "yes","no")) -> new_sites

ggplot(new_sites %>% filter (gen%%1== 0 )) + 
  geom_histogram(aes(gen),stat="count",fill="black") +
  facet_wrap(~location) +
  theme_bw() + 
  ylab("GBS Sample size") + 
  xlab("Generation")


###Only LBJ
new_sites <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.BC1.LBJ.missing50.dp5.10mac.bi.newalleles.txt",header=T,
                       colClasses = c("factor","factor","numeric","numeric","numeric","numeric"))
#new_sites <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.BC1.LBJ.missing50.dp5.3mac.bi.newalleles.txt",header=T,
#                        colClasses = c("factor","factor","numeric","numeric","numeric","numeric","numeric"))

new_sites  %>%
  mutate(gen = gsub("2.5","3",gen)) %>%
  mutate(gen = gsub("7.5","8",gen)) %>%
  filter(location == "LBJ" | gen == "1") %>%
  ggplot(.) + 
  geom_boxplot(aes(x=as.factor(gen),y=percent_new*100),color="lightgrey", outlier.size=0) +
  geom_jitter(aes(x=as.factor(gen),y=percent_new*100),width=0.2) + 
  theme_few() +
  ylab("Percent of variants not found in BC1 generation") + 
  xlab("Generation") +
  labs(tag = "B") +
  scale_x_discrete(limits = c("1", "2", "3", "4","5","6","7","8"))-> plot.2
  
gen_1_contam <- tibble(gen=as.character(1),contaminated=as.character("no"),sum_type=as.numeric(100))

new_sites %>% filter(gen != "1") %>%
  mutate(gen = gsub("2.5","3",gen)) %>%
  mutate(gen = gsub("7.5","8",gen)) %>%
  filter(location == "LBJ") %>%
  mutate(contaminated = ifelse(percent_new > 0.005, "yes","no")) %>% 
  group_by(gen, contaminated) %>%
  summarize(sum_type = n()) %>%  
  ungroup() %>%
  rbind(., gen_1_contam) %>% 
  ggplot(aes(x=gen,y=sum_type,fill=contaminated))+ #define the x and y 
  geom_bar(stat = "identity",position="fill") + #make the stacked bars
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1",name="Admixture",labels=c("No","Yes")) +
  ylab("Percent") + xlab("Generation") +theme_few() +
  labs(tag = "A") +
  theme(legend.position="bottom") +
  scale_x_discrete(limits = c("1", "2", "3", "4","5","6","7","8")) -> plot.1
  
gen_1_contam_sd <- tibble(gen=as.character(1),contam_sd=as.numeric(0))


new_sites %>% filter(gen != "1") %>%
  mutate(gen = gsub("2.5","3",gen)) %>%
  mutate(gen = gsub("7.5","8",gen)) %>%
  filter(location == "LBJ") %>%
  group_by(gen) %>%
  summarize(contam_sd = sd(percent_new*100)) %>% 
  rbind(.,gen_1_contam_sd) %>%
  ggplot(.,aes(x=as.factor(gen),y=contam_sd)) +
  geom_bar(stat="identity") +
  xlab("Generation") +theme_few() +
  ylab("Admixture standard deviation") +
  labs(tag = "C") +
  scale_x_discrete(limits = c("1", "2", "3", "4","5","6","7","8")) +
  geom_segment(aes(x=0.5,xend=1.5,y=0,yend=0)) -> plot.3

pdf("texanus.gatk.GBS2.BC1.LBJ.missing50.dp5.10mac.bi.contamination.pdf",height=12,width=6)
#pdf("texanus.gatk.GBS2.BC1.LBJ.missing50.dp5.3mac.bi.contamination.pdf",height=12,width=6)

grid.arrange(plot.1, plot.2,plot.3, ncol = 1)
dev.off()

new_sites %>% filter(gen != "1") %>% summarize(max_count = max(total_sites),min_count = min(total_sites),mean_count = mean(total_sites))
