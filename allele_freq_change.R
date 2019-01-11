library(stringr)
library(tidyverse)
library(broom)
library(RColorBrewer)
colors <-brewer.pal(n = 8, name = "Set1")
poplist <- c("BFL","HCC","KPC","LBJ")
folder <- "/home/owens/working/texanus/allele_freq"
files <- list.files(folder, pattern='.frq')
data <- data.frame(chr=factor(),
                   pos=numeric(),
                   freq=numeric(),
                   gen=numeric(),
                   location=factor())
for (i in 1:length(files)){
  tmp <- read.delim(paste(folder,files[i],sep="/"),header=F,skip=1)
  colnames(tmp) <- c("chr","pos","n_alleles","n_ind","ref","freq")
  tmp %>% select(chr, pos, freq) -> tmp
  tmp2 <- gsub("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.","",files[i],)
  tmp2 <- gsub("_list.frq","",tmp2)
  info <- str_split(tmp2[[1]][[1]],"_")
  location <- info[[1]][1]
  gen <- as.numeric(info[[1]][[2]])
  tmp$gen <- gen
  tmp$location <- location
  data <- rbind(data, tmp)
}
gen1 <- read.delim(paste(folder,"texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.all_1_list.frq",sep="/"),header=F,skip=1)
colnames(gen1) <- c("chr","pos","n_alleles","n_ind","ref","freq")
gen1 %>% select(chr, pos, freq) -> gen1
for (i in c("HCC","LBJ","BFL","KPC")){
tmp <- gen1
tmp$gen <- as.numeric(1)
tmp$location <- i
data <- rbind(data, tmp)
}
gen1 %>% select(chr, pos, freq) -> gen1
colnames(gen1) <- c("chr","pos","gen1_freq")
data <- merge(data, gen1)
data$change <- data$freq - data$gen1_freq

data <- merge(data, marker_identity)
chrlist <- data$chr %>% unique()
#wgs_snps <- read.delim("/home/owens/working/texanus/WGS/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.15-35af.WGScompared.txt",header=T,comment.char="#")
#wgs_snps %>% select(chr, pos) -> wgs_snps_pos
#inner_join(data, wgs_snps_pos) %>% filter(gen1_freq > 0.15, gen1_freq < 0.35) -> data_wgsfiltered

#tmp <- merge(data, lg,all.x =TRUE)


#ld_snps <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.multipopld.txt",header=T)
#ld_snps %>%
#  count(chr, haplotype) %>% arrange(desc(n)) -> ld_counts
#ld_snps <- merge(ld_snps, ld_counts)

#data <- merge(data, ld_snps)

###
pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.deballeles.change.pdf",height=20,width=20)
for (i in 1:4){
  print(
  data %>% filter(parent == "Debilis") %>%
    filter(location == poplist[i]) %>% 
    # filter(chr == "HanXRQChr14") %>%
    ggplot(.,aes(x=pos,y=freq)) + 
    geom_point() +
    facet_grid(chr~gen) +
    ggtitle(poplist[i])
  )
}
for (i in 1:17){
  print(
    data %>% filter(parent == "Debilis") %>%
      filter(location != "all") %>% 
      filter(chr == chrlist[i]) %>%
      ggplot(.,aes(x=pos,y=freq)) + 
      geom_point() +
      facet_grid(location~gen) +
      ggtitle(chrlist[i])
  )
}

dev.off()






data %>% filter(chr !="HanXRQChr00c0494", chr !="HanXRQChr00c0041") -> data
data %>% filter(gen1_freq > 0.15, gen1_freq < 0.35) %>%
  filter(location == "LBJ") %>% 
  filter(chr == "HanXRQChr01") %>%
  ggplot(., aes(x=gen, y=change))  + geom_point() + geom_smooth(method="lm")

pdf("Texanus_allele_freq_change_haplotype_example.v0.pdf",width=6,height=6)
data %>% filter(parent == "Debilis") %>%
  filter(location == "BFL") %>% 
  filter(gen == "14") %>%
  filter(chr == "HanXRQChr17") %>%
 # filter(haplotype == 2) %>%
  ggplot(.,aes(x=pos,y=change)) + 
  geom_point() 
  #geom_smooth()
dev.off()

data %>% filter(parent == "Debilis") %>%
  filter(location != "all") %>%
  #filter(location == "BFL") %>%
  #filter(chr == "HanXRQChr02") %>%
  #filter(haplotype == 2) %>%
  #filter(chr == "HanXRQChr06") %>%
  #filter(haplotype == 274) %>%
  #filter(chr == "HanXRQChr17") %>%
  ggplot(.,aes(x=gen,y=change,group=gen)) + geom_boxplot() +facet_grid(chr~location)

pdf("Texanus_allele_freq_change_v0.6.pdf",width=12,height=8)
data %>% filter(parent == "Debilis") %>%
  filter(location == "BFL",gen ==13) %>% 
  ggplot(.) + geom_point(aes(x=pos,y=change),color=colors[1],alpha=0.2) + 
  geom_point(data=data_wgsfiltered %>% filter(gen1_freq > 0.2, gen1_freq < 0.3,location=="HCC", gen==9),aes(x=pos,y=change),color=colors[2],alpha=0.2)+
  geom_point(data=data_wgsfiltered %>% filter(gen1_freq > 0.2, gen1_freq < 0.3,location=="KPC", gen==8),aes(x=pos,y=change),color=colors[3],alpha=0.2)+
  geom_point(data=data_wgsfiltered %>% filter(gen1_freq > 0.2, gen1_freq < 0.3,location=="LBJ", gen==7.5),aes(x=pos,y=change),color=colors[4],alpha=0.2)+
  geom_smooth(aes(x=pos,y=change),color=colors[1]) +
  geom_smooth(data=data_wgsfiltered %>% filter(gen1_freq > 0.2, gen1_freq < 0.3,location=="HCC", gen==9),aes(x=pos,y=change),color=colors[2])+
  geom_smooth(data=data_wgsfiltered %>% filter(gen1_freq > 0.2, gen1_freq < 0.3,location=="KPC", gen==8),aes(x=pos,y=change),color=colors[3])+
  geom_smooth(data=data_wgsfiltered %>% filter(gen1_freq > 0.2, gen1_freq < 0.3,location=="LBJ", gen==7.5),aes(x=pos,y=change),color=colors[4])+
  theme_bw() +  
  geom_hline(yintercept=0) + 
  facet_wrap(~chr) +
  ggtitle("Change in allele frequency from generation 1 to last generation")+ 
  theme(axis.text.x=element_blank()) 
dev.off()


data %>% group_by(chromo,position) %>%
  filter(chromo == "HanXRQChr01") %>%
  mutate(gen1 = knownEM[gen == "1"]) %>% head()

data %>% filter(gen_1_freq > 0.2, gen_1_freq < 0.3) %>%
  filter(gen >6) %>% 
  filter(change > 0.1) %>% 
  ggplot(.) + geom_point(aes(x=position,y=change,color=location)) + facet_wrap(~chromo)


regressions <- data %>% filter(location != "all") %>%
  filter(gen_1_freq > 0.2, gen_1_freq < 0.3) %>%
  group_by(location, chromo, position) %>%
  do(fit = lm(change ~ 0 + gen, .))

regression.fit <- regressions %>% tidy(fit) 

regression.fit %>% filter(term == "gen") %>% filter(statistic != "NaN") %>%
  #filter(location == "BFL") %>%
  ggplot(aes(x=position,y=estimate,color=location)) + 
  #geom_point() + 
  geom_smooth(span = 1000) +
  facet_wrap(~chromo) +
  geom_hline(yintercept=0)

data_wgsfiltered %>% filter(gen != 1, chr == "HanXRQChr02") %>% 
  ggplot(.,aes(change, colour = as.factor(gen))) + 
  geom_freqpoly() + facet_wrap(~location)

data_wgsfiltered %>% filter(gen != 1, chr == "HanXRQChr02", pos <100000000 ) %>% 
  ggplot(.,aes(y=freq, x= gen)) + 
  geom_jitter() + facet_wrap(~location) + geom_smooth(method=lm)




  