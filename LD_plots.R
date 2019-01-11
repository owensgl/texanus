gen_1 <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.all_1.geno.ld",header=T)
gen_14 <- read.delim("/home/owens/working/texanus/texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi.BFL_14.geno.ld",header=T)

gen_1$dist <- abs(gen_1$POS1 - gen_1$POS2)
gen_14$dist <- abs(gen_14$POS1 - gen_14$POS2)

gen_14 %>% head(100000) %>%
  filter(R.2 != "-nan") %>%
  ggplot(aes(x=dist,y=R.2)) + geom_point() + geom_smooth()


colnames(gen_1) <- c("CHR","POS1","POS2","N_INDV","GEN1R2","dist")
colnames(gen_14) <- c("CHR","POS1","POS2","N_INDV","GEN14R2","dist")

both <- merge(gen_1,gen_14)

both %>% head(10000) %>%
  ggplot(., aes(x=GEN1R2,y=GEN14R2,color=dist)) + geom_point() + geom_smooth()
