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

k = 1
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
    rel_freq = n_alleles/mean(total_local),
    total_alleles = mean(total_local)
  ) %>%  
  group_by(location, gen, chr,start) %>%
  rbind(., gen1) %>%
  mutate(n_parents = n()) %>% 
  filter(n_parents == 4) -> window_frequency

window_frequency %>% filter(chr == "HanXRQChr01", 
                            start == 0, 
                            location == "HCC",
                            gen == "6") %>%
  ungroup() %>%
  dplyr::select(rel_freq) %>% .$rel_freq


