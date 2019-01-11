#plot FST between texanus and central annuus
library(tidyverse)
library(ggthemes)
library(ggrastr)

directory <- "/home/owens/working/texanus/WGS"
# 
# ANNCEN_ANNTEX <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.ANNCEN_ANNTEX.fst.txt",sep=""))
# ANNCEN_ANNTEX$comparison <- "ANNCEN_ANNTEX"
# ANNCEN_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.ANNCEN_DEB.fst.txt",sep=""))
# ANNCEN_DEB$comparison <- "ANNCEN_DEB"
# ANNTEX_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.ANNTEX_DEB.fst.txt",sep=""))
# ANNTEX_DEB$comparison <- "ANNTEX_DEB"
# ANNTEX_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.ANNTEX_ARG.fst.txt",sep=""))
# ANNTEX_ARG$comparison <- "ANNTEX_ARG"
# ANNCEN_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.ANNCEN_ARG.fst.txt",sep=""))
# ANNCEN_ARG$comparison <- "ANNCEN_ARG"

ANNCEN_ANNTEX <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412.ANNCEN_ANNTEX.fst.txt.gz",sep=""))
ANNCEN_ANNTEX$comparison <- "ANNCEN_ANNTEX"
ANNCEN_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412.ANNCEN_DEB.fst.txt.gz",sep=""))
ANNCEN_DEB$comparison <- "ANNCEN_DEB"
ANNTEX_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412.ANNTEX_DEB.fst.txt.gz",sep=""))
ANNTEX_DEB$comparison <- "ANNTEX_DEB"
ANNTEX_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412.ANNTEX_ARG.fst.txt.gz",sep=""))
ANNTEX_ARG$comparison <- "ANNTEX_ARG"
ANNCEN_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412.ANNCEN_ARG.fst.txt.gz",sep=""))
ANNCEN_ARG$comparison <- "ANNCEN_ARG"

# fst <- rbind(ANNCEN_ANNTEX, ANNCEN_DEB, ANNTEX_DEB, ANNCEN_ARG, ANNTEX_ARG)
fst <- rbind(ANNCEN_ANNTEX, ANNCEN_DEB, ANNTEX_DEB,ANNCEN_ARG, ANNTEX_ARG)
window_size = 1000000
pdf("Texanus.tranche90.snp.gwas.90.bi.multiFst.v2.HA412.pdf",width=12,height=8)
fst %>%
  filter(comparison == "ANNCEN_ANNTEX") %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>%
  filter(count > 10) %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window/window_size,y=window_fst)) + 
  geom_point_rast(data=fst %>% filter(comparison == "ANNCEN_ANNTEX"), aes(x=pos/window_size,y=Fst),alpha=0.1)+
  facet_wrap(~chr) +
  theme_few() + geom_line(aes(),alpha=1,color="red") + 
  xlab("MB") + ggtitle("Annuus Fst: Central <-> Texas") 


fst %>%
  filter(comparison != "ANNCEN_DEB", comparison != "ANNTEX_DEB") %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom), count=n()) %>%
  filter(count > 10) %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window/window_size,y=window_fst,group=comparison)) + facet_wrap(~chr) +
  theme_few() + geom_line(aes(color=comparison),alpha=0.7) + 
  scale_color_brewer(palette = "Set1",name="Comparison") + xlab("MB") +
  ggtitle("Annuus - Argophyllus Fst")

fst %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>% 
  filter(count > 10) %>%
  select(-count) %>%
  spread(comparison, window_fst) %>% 
  mutate(Fstdif = (ANNCEN_ARG - ANNTEX_ARG)) %>%
  filter(!is.na(ANNCEN_ANNTEX)) %>%
  filter(!is.na(Fstdif)) %>% 
  ggplot() + facet_wrap(~chr) +
  theme_few() + geom_line(aes(x=window/window_size,y=ANNCEN_ANNTEX),color="red") + 
  geom_line(aes(x=window/window_size,y=Fstdif),color="black") +
  geom_hline(yintercept=0,linetype="dotted",color="grey") + xlab("MB") +
  ylab("Fst or Fst difference") +
  ggtitle("ANNCEN_ANNTEX Fst = red\nANNCEN_ARG - ANNTEX_ARG Fst = black")


fst %>%
  filter(comparison != "ANNCEN_ARG", comparison != "ANNTEX_ARG") %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>%
  filter(count > 10) %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window/window_size,y=window_fst,group=comparison)) + facet_wrap(~chr) +
  theme_few() + geom_line(aes(color=comparison),alpha=0.7) + 
  scale_color_brewer(palette = "Set1",name="Comparison") + xlab("MB") +
  ggtitle("Annuus - Debilis Fst")

fst %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>% 
  filter(count > 10) %>%
  select(-count) %>%
  spread(comparison, window_fst) %>% mutate(Fstdif = ANNCEN_DEB - ANNTEX_DEB) %>% 
  filter(!is.na(ANNCEN_ANNTEX)) %>%
  filter(!is.na(Fstdif)) %>%
  ggplot() + facet_wrap(~chr) +
  theme_few() + geom_line(aes(x=window/window_size,y=ANNCEN_ANNTEX),color="red") + 
  geom_line(aes(x=window/window_size,y=Fstdif),color="black") +
  geom_hline(yintercept=0,linetype="dotted",color="grey") + xlab("MB") + 
  ylab("Fst or Fst difference") +
  ggtitle("ANNCEN_ANNTEX Fst = red\nANNCEN_DEB - ANNTEX_DEB Fst = black")
dev.off()

window_size = 100000
pdf("Texanus.tranche90.snp.gwas.90.bi.multiFst.debwindows.v1.pdf",height=5,width=9)
fst %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom)) %>% 
  spread(comparison, window_fst) %>% mutate(Fstdif = ANNCEN_DEB - ANNTEX_DEB) %>%
  filter(!is.na(ANNCEN_ANNTEX)) %>%
  filter(!is.na(Fstdif)) %>%
  mutate(outlier = case_when(ANNCEN_ANNTEX > 0.4 & Fstdif > 0.4 ~ "outlier",
                             TRUE ~ "non-outlier")) %>%
  ggplot(aes(x=ANNCEN_ANNTEX, y=Fstdif)) + geom_point(aes(color=outlier)) +
  theme_few() + geom_smooth(method="lm") +
  scale_color_manual(values=c("black","red")) +
  geom_hline(yintercept=0,color="grey") +
  xlab("ANNCEN_ANNTEX Fst") +
  ylab("ANNCEN_DEB - ANNTEX_DEB Fst") +
  ggtitle(paste(window_size,"windows"))
dev.off()

fst %>%
  mutate(window = floor(pos/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom)) %>% 
  spread(comparison, window_fst) %>% mutate(Fstdif = ANNCEN_DEB - ANNTEX_DEB) %>%
  filter(!is.na(ANNCEN_ANNTEX)) %>%
  filter(!is.na(Fstdif)) %>%
  mutate(outlier = case_when(ANNCEN_ANNTEX > 0.4 & Fstdif > 0.4 ~ "outlier",
                             TRUE ~ "non-outlier")) %>%
  filter(outlier == "outlier") %>%
  mutate(window_end = window + window_size) %>%
  select(chr, window, window_end, everything()) %>%
  write_tsv("Texanus.tranche90.snp.gwas.90.bi.fst.doubleoutlier.v1.txt")


