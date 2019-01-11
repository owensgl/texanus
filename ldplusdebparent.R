library(stringr)
library(tidyverse)
library(broom)
library(RColorBrewer)

folder <- "/home/owens/working/texanus"
ld.test <- read.delim(paste(folder,"texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.all_1.debcompared.txt",sep="/"),header=T)

ld.test %>%
ggplot(.) + geom_histogram(aes(x=percent,fill=deb_identified)) +
  facet_wrap(~chr) + geom_vline(xintercept=0.25)

ld.test %>% filter(percent > 0.25) -> deb.alleles

ld.test$col <- cut(ld.test$percent,
               breaks = c(0, 0.1, Inf),
               labels = c("low", "high"))
ld.test %>%
ggplot(.) + geom_jitter(aes(x=pos,y=1,col=col)) + facet_wrap(~chr)
