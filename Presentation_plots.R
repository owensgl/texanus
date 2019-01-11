#Presentation plots

new_sites <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.newalleles.txt",header=T,
                        colClasses = c("factor","factor","numeric","numeric","numeric","numeric"))
pdf("Texanus_sample_size.pdf")
ggplot(new_sites %>% filter (gen%%1== 0 )) + 
  geom_histogram(aes(gen),stat="count",fill="black") +
  facet_wrap(~location) +
  theme_bw() + 
  ylab("GBS Sample size") + 
  xlab("Generation")
dev.off()

gen1 %>% filter(location == "BFL") 
  ggplot(.) + geom_point(aes(x=)) 
