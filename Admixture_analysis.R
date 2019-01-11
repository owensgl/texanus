#Loading all individual faststructure analyses
library(tidyverse)

directory <- "/home/owens/working/texanus/WGS"

pdf("Texanus.tranche90.snp.gwas.90.bi.0p2ld.admix.thin.admixture.pdf",height=15,width=15)

Q <- read_table(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.0p2ld.admix.thin.6.Q",sep=""),col_names=c("Q1","Q2","Q3","Q4","Q5","Q6"))
samples <- read_table(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.admix.samplelist.txt",sep=""),col_names="name")
all_info <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv")
subspecies_info <- read_tsv(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.admix.sampleinfo.txt",sep=""),col_names=c("name","subspecies"))

data <- cbind(samples,Q)

inner_join(data, all_info) %>% 
  inner_join(., subspecies_info) %>% 
  gather(Q, value, Q1:Q6) %>% mutate(subspecies = gsub("-","ANN_TEX",subspecies)) %>%
  ggplot(aes(x = name, y = value, fill = factor(Q)))+
  geom_bar(stat = "identity", width = 1.1) +facet_wrap(~subspecies,scale="free_x",ncol=1) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_brewer(palette = "Set1",name="Group") +
  ggtitle("Admixture results")
    
dev.off()

