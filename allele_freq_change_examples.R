pdf("Texanus_example_freq_change_haplotypes.v0.pdf")
data %>% filter(gen1_freq > 0.15, gen1_freq < 0.35) %>%
  filter(location != "all") %>%
  filter(chr == "HanXRQChr02") %>%
  filter(haplotype == 2) %>%
  ggplot(.,aes(x=gen,y=change,group=gen)) + geom_boxplot() +facet_wrap(~location) +
  ggtitle("Chr2")

data %>% filter(gen1_freq > 0.15, gen1_freq < 0.35) %>%
  filter(location != "all") %>%
  filter(chr == "HanXRQChr06") %>%
  filter(haplotype == 274) %>%
  ggplot(.,aes(x=gen,y=change,group=gen)) + geom_boxplot() +facet_wrap(~location) +
  ggtitle("Chr6")

data %>% filter(gen1_freq > 0.15, gen1_freq < 0.35) %>%
  filter(location != "all") %>%
  filter(chr == "HanXRQChr12") %>%
  filter(haplotype == 191) %>%
  ggplot(.,aes(x=gen,y=change,group=gen)) + geom_boxplot() +facet_wrap(~location) +
  ggtitle("Chr12")

data %>% filter(gen1_freq > 0.15, gen1_freq < 0.35) %>%
  filter(location != "all") %>%
  filter(chr == "HanXRQChr15") %>%
  filter(haplotype == 522) %>%
  ggplot(.,aes(x=gen,y=change,group=gen)) + geom_boxplot() +facet_wrap(~location) +
  ggtitle("Chr15")

dev.off()