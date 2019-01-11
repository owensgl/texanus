#Plotting Fst 
library(tidyverse)
fst.wind <- read.delim("/home/owens/working/texanus/WGS/texanusWGS.gatk.MAF5.min2.500kbSW.fst.txt",header=T)

fst.wind$Chr <- gsub("HanXRQChr","",fst.wind$Chr)
fst.wind <- fst.wind %>% filter(., !grepl("00",Chr)) 

fst.wind %>% 
  ggplot(.,aes(x=MidPos,y=Fst)) + geom_line()+
  facet_wrap(~Chr)

fst <- read.delim("/home/owens/working/texanus/WGS/texanusWGS.gatk.MAF5.min2.fst.txt",header=T)

fst %>%
  ggplot(.,aes(Fst)) + geom_histogram()
