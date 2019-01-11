library(tidyverse)
library(viridis)
library(broom)
library(segmented)
chrlist <- read.delim("/home/owens/working/texanus/rqtl/chrlist.txt",header=F)

#window <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.10mb.txt",header=T,
#                        colClasses = c("factor","factor","factor","numeric","numeric","factor","numeric")))
window <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.txt",header=T,
                     colClasses = c("character","character","character","numeric","numeric","character","character"))


info <- read.delim("/home/owens/working/texanus/texanus.sampleinfo.txt",header=T)
info %>% dplyr::select(ID, gen, location) -> info
colnames(info) <- c("sample","gen","location")
poplist <- c("BFL","LBJ","HCC","KPC")
window <- merge(window, info)
window$chr <- as.factor(window$chr)
window %>% filter(call != "Draw",gen == 1) %>%
  group_by(gen,chr,start) %>%
  mutate(total_local = sum(as.numeric(as.character(call)))) %>%
  
  group_by(gen, parent, chr, start, end) %>%
  summarize(n_genotypes = n() * 2, 
            n_alleles = sum(as.numeric(as.character(call))),
            freq = n_alleles / n_genotypes,
            rel_freq = n_alleles/mean(total_local)) -> gen1
gen1.tmp <- gen1
gen1$location <- "BFL"
for (i in 2:4){
  tmp <- gen1.tmp
  tmp$location <- poplist[i]
  gen1 <- rbind(gen1,tmp)
}
pdf(
  "texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.freqpergen.pdf"
)
for (k in 1:nrow(chrlist)) {
  print(
    window %>% filter(call != "Draw") %>%
      filter(gen != 1) %>%
      filter(gen != 7.5) %>%
      group_by(location,gen,chr,start) %>%
      mutate(total_local = sum(as.numeric(as.character(call)))) %>%
      group_by(location, gen, parent, chr, start) %>%
      summarize(
        n_genotypes = n() * 2,
        n_alleles = sum(as.numeric(as.character(call))),
        freq = n_alleles / n_genotypes,
        rel_freq = n_alleles/mean(total_local)
      ) %>%
      rbind(., gen1) %>%
      filter(chr == as.character(chrlist[k, ])) %>%
      ggplot(.) +
      geom_line(aes(
        x = gen,
        y = freq,
        group = interaction(start, parent),
        color = parent
      )) +
      facet_grid(location ~ start) +
      theme_minimal() +
      scale_color_brewer(palette = "Set1") +
      ylab("Haplotype Frequency") +
      xlab("Generation") +
      ggtitle(paste(
        "20MB Windowed Haplotype frequency", chrlist[k, ], sep = " "
      ))
  )
}
dev.off()

pdf(
  "texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.relfreqpergen.pdf"
)
for (k in 1:nrow(chrlist)) {
  print(
    window %>% filter(call != "Draw") %>%
      filter(gen != 1) %>%
      filter(gen != 7.5) %>%
      group_by(location,gen,chr,start) %>%
      mutate(total_local = sum(as.numeric(as.character(call)))) %>% 
      group_by(location, gen, parent, chr, start) %>%
      summarize(
        n_genotypes = n() * 2,
        n_alleles = sum(as.numeric(as.character(call))),
        freq = n_alleles / n_genotypes,
        rel_freq = n_alleles/mean(total_local)
      ) %>%
      group_by(location,gen, chr,start) %>%
      rbind(., gen1) %>%
      mutate(n_parents = n()) %>%
      filter(n_parents == 4) %>%
      filter(chr == as.character(chrlist[k, ])) %>%
      ggplot(.) +
      geom_line(aes(
        x = gen,
        y = rel_freq,
        group = interaction(start, parent),
        color = parent
      )) +
      facet_grid(location ~ start) +
      theme_minimal() +
      scale_color_brewer(palette = "Set1") +
      ylab("Relative Haplotype Frequency") +
      xlab("Generation") +
      ggtitle(paste(
        "20MB Windowed Haplotype frequency", chrlist[k, ], sep = " "
      ))
  )
}
dev.off()




######Linear model for haplotype frequency
window %>% filter(call != "Draw") %>%
  filter(gen != 1) %>%
  filter(gen != 7.5) %>%
  group_by(location,gen,chr,start) %>%
  mutate(total_local = sum(as.numeric(as.character(call)))) %>%
  group_by(location, gen, parent, chr, start,end) %>%
  summarize(
    n_genotypes = n() * 2,
    n_alleles = sum(as.numeric(as.character(call))),
    freq = n_alleles / n_genotypes,
    rel_freq = n_alleles/mean(total_local)
  ) %>%
  rbind(., gen1) %>%
  ungroup() %>%
  group_by(location, parent, chr, start,end) %>%
  do(tidy(lm(I(freq - 0.25) ~ 0 + gen, .))) -> window.lm


window.lm$chr <- as.factor(window.lm$chr)
levels(window.lm$chr) <- 1:17

pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.freqpergenslopes.pdf",
    width=15,height=2.5)
window.lm %>% 
  filter(term == "gen") %>% 
  ggplot(.) + 
  geom_segment(aes(x=start, xend=end,y=as.factor(parent),yend=as.factor(parent),color=estimate),size=2) +
  scale_color_distiller(palette = "RdYlBu", limits=c(-.08,.08),name="Slope") +
  facet_grid(location~chr) + theme_minimal() +
  ylab("Location") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Change in haplotype frequency in 20MB windows")

window.lm %>% 
  filter(term == "gen") %>% 
  filter(parent == "deb") %>%
  ggplot(.) + 
  geom_segment(aes(x=start, xend=end,y=as.factor(location),yend=as.factor(location),color=estimate),size=2) +
  scale_color_distiller(palette = "RdYlBu", limits=c(-.05,.05),name="Slope") +
  facet_grid(~chr) + theme_minimal() +
  ylab("Location") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Change in debilis haplotype frequency in 20MB windows")
dev.off()
  
  
#### Relative frequency linear models
window %>% filter(call != "Draw") %>%
  filter(gen != 1) %>%
  filter(gen != 7.5) %>%
  group_by(location,gen,chr,start) %>%
  mutate(total_local = sum(as.numeric(as.character(call)))) %>%
  group_by(location, gen, parent, chr, start,end) %>%
  summarize(
    n_genotypes = n() * 2,
    n_alleles = sum(as.numeric(as.character(call))),
    freq = n_alleles / n_genotypes,
    rel_freq = n_alleles/mean(total_local)
  ) %>%
  rbind(., gen1) %>%
  ### For relative freq calculations
  group_by(location,gen, chr,start) %>%
  mutate(n_parents = n()) %>% 
  filter(n_parents == 4) %>% 
  ###
  ungroup() %>%
  group_by(location, parent, chr, start,end) %>%
  do(tidy(lm(I(rel_freq - 0.25) ~ 0 + gen, .))) -> windowrel.lm


windowrel.lm$chr <- as.factor(windowrel.lm$chr)
levels(windowrel.lm$chr) <- 1:17

pdf("texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.parentagewindows.20mb.relfreqpergenslopes.pdf",
    width=15,height=2.5)
windowrel.lm %>% 
  filter(term == "gen") %>% 
  ggplot(.) + 
  geom_segment(aes(x=start, xend=end,y=as.factor(parent),yend=as.factor(parent),color=estimate),size=2) +
  scale_color_distiller(palette = "RdYlBu", limits=c(-.081,.081),name="Slope") +
  facet_grid(location~chr) + theme_minimal() +
  ylab("Location") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Change in haplotype relative frequency in 20MB windows")

windowrel.lm %>% 
  filter(term == "gen") %>% 
  filter(parent == "deb") %>% 
  ggplot(.) + 
  geom_segment(aes(x=start, xend=end,y=as.factor(location),yend=as.factor(location),color=estimate),size=2) +
  scale_color_distiller(palette = "RdYlBu", limits=c(-.054,.054),name="Slope") +
  facet_grid(~chr) + theme_minimal() +
  ylab("Location") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Change in debilis haplotype relative frequency in 20MB windows") +
  theme(panel.spacing.y=unit(0.1, "lines"))
  
dev.off()





####Trying to do segmented linear models
window %>% filter(call != "Draw") %>%
  filter(gen != "1.0") %>%
  filter(gen != 7.5) %>%
  group_by(location, gen, parent, chr, start,end) %>%
  summarize(
    n_genotypes = n() * 2,
    n_alleles = sum(as.numeric(as.character(call))),
    freq = n_alleles / n_genotypes
  ) %>%
  rbind(., gen1) %>%
  ungroup() %>%
  group_by(location, parent, chr, start,end) %>%
  filter(chr == "HanXRQChr04") %>%
  do(tidy(segmented(lm(I(freq - 0.25) ~ 0 + gen, .),seg.Z =~gen))) -> tmp

 
 segmented.mod <- segmented(tmp.lm,seg.Z =~ gen)
 plot(tmp$gen, tmp$freq, ylim=c(-1,1))
 plot(segmented.mod, add=T)

 
abba <- read.delim("texanus.abba.centralann.stex.results",header=F,sep=" ")
colnames(abba)  <- c("chr","start","close_end","denom","sites", "Dstat")
merge(window.lm,abba) -> window.lm.abba

pdf("texanus_Dstat_prelim.v1.pdf")
window.lm.abba %>% filter(parent == "deb") %>%
ggplot(.,aes(x=Dstat,y=estimate,color=location)) + geom_point() + 
  geom_smooth(method="lm") + ylab("experimental allele freq slope")

window.lm.abba %>% filter(parent == "deb") %>%
  ggplot(.,aes(x=start,y=Dstat)) + geom_point(aes()) + 
  facet_wrap(~chr)

window.lm.abba %>% filter(parent == "deb") %>%
  ggplot(.,aes(Dstat)) + geom_histogram()
dev.off()
       
       


window %>% filter(call != "Draw") %>%
  filter(gen != 1) %>%
  filter(gen != 7.5) %>%
  filter(location == "LBJ") %>%
  group_by(gen,chr,start) %>%
  mutate(total_local = sum(as.numeric(as.character(call)))) %>%
  group_by(location, gen, parent, chr, start) %>%
  summarize(
    n_genotypes = n() * 2,
    n_alleles = sum(as.numeric(as.character(call))),
    freq = n_alleles / n_genotypes,
    rel_freq = n_alleles/mean(total_local)
  ) %>%
  rbind(., gen1) %>%
  filter(parent == "deb") %>%
  ggplot(.) +
  geom_line(aes(
    x = (start+10000000)/1000000,
    y = freq,
    group = gen,
    color = gen
  )) +
  facet_wrap(~chr,nrow=3,scales="free_x") +
  theme_few() +
  scale_color_distiller(palette = "Spectral") +
  ylab("Debilis Frequency") +
  xlab("MB")

window %>% filter(call != "Draw") %>%
  filter(gen != 1) %>%
  filter(location == "LBJ") %>% 
  filter(parent == "deb") %>%
  group_by(sample,gen, call) %>%
  summarize(count=n()) %>%
  spread(call,count,fill=0) %>% 
  mutate(percent_deb = ((`2`* 2)+`1`) / ((`0` + `1` + `2`) * 2)) %>% 
  ggplot(.,aes(x=as.factor(gen),y=percent_deb)) + geom_boxplot() +
  theme_few() + xlab("Generation") + ylab("Proportion debilis counted by 20MB window") 


window %>% filter(call != "Draw") %>%
  filter(gen != 1) %>%
  filter(location == "LBJ") %>% 
  filter(parent == "deb") %>%
  group_by(sample,gen, call) %>%
  summarize(count=n()) %>%
  spread(call,count,fill=0) %>% 
  mutate(percent_deb = ((`2`* 2)+ `1`) / ((`0` + `1` + `2`) * 2)) %>% 
  filter(gen == 7) %>%View()
  mutate(any_debilis = case_when(percent_deb >  ~ "Yes",
                                 TRUE ~ "No")) %>% 
  ggplot(.,aes(x=as.factor(gen),y=any_debilis,fill=any_debilis)) + geom_bar(stat="identity") +
  theme_few() + xlab("Generation") + ylab("Percent debilis by windows")
 