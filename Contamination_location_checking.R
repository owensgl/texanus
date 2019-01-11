#This script is for plotting the contamination allele frequency.


directory <- "/home/owens/working/texanus"

contam <- read_tsv(paste(directory,"/texanus.gatk.LBJ.contamination.all.txt",sep=""))

pdf("texanus.gatk.LBJ.contamination.all.pdf")
contam %>% 
  filter(!is.na(ANN_42)) %>% filter(!is.na(ANN_43)) %>% 
  filter(ANN_42 != ANN_43) %>% 
  mutate(unique = case_when(ANN_42 == 0 ~ "ANN_43",
                            ANN_43 == 0 ~ "ANN_42",
                            TRUE ~ "neither")) %>% 
  filter(unique != "neither") %>%
  ggplot(aes(unique)) + geom_bar(stat="count") +
  theme_bw() + ylab("Count") + xlab("Contaminated alleles in one population")
dev.off()
total_tests <- contam %>% 
  filter(!is.na(ANN_42)) %>% filter(!is.na(ANN_43)) %>% 
  filter(ANN_42 != ANN_43) %>% 
  mutate(unique = case_when(ANN_42 == 0 ~ "ANN_43",
                            ANN_43 == 0 ~ "ANN_42",
                            TRUE ~ "neither")) %>%
  filter(unique != "neither") %>% nrow()

total_ANN_42 <- contam %>% 
  filter(!is.na(ANN_42)) %>% filter(!is.na(ANN_43)) %>% 
  filter(ANN_42 != ANN_43) %>% 
  mutate(unique = case_when(ANN_42 == 0 ~ "ANN_43",
                            ANN_43 == 0 ~ "ANN_42",
                            TRUE ~ "neither")) %>%
  filter(unique == "ANN_42") %>% nrow()

binom.test( total_ANN_42,
            total_tests,
           p=0.5,
           alternative="two.sided")



contam %>% 
  filter(!is.na(ANN_42)) %>% filter(!is.na(ANN_43)) %>% 
  filter(ANN_42 != ANN_43) %>% 
  mutate(unique = case_when(ANN_42 == 0 ~ "ANN_43",
                            ANN_43 == 0 ~ "ANN_42",
                            TRUE ~ "neither")) %>%
  filter(unique == "neither") %>%
  mutate(diff = ANN_42 - ANN_43) %>%
  summarize(mean_diff = mean(diff))

