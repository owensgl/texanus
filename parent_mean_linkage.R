
basename <- "texanus.gatk.GBS2.noANN.missing50.dp5.10mac.bi"
parents <- c("deb","ann1","ann2","ann3")
parental_linkage <- data.frame(chr=as.character(),pos=as.numeric(),n_sites=as,numeric(),
                               mean_r2=as.numeric(),location=as.character(),gen=as.numeric(),
                               parents=as.character())
for (i in 1:4){
  genotypes %>% filter(location == poplist[i]) %>% select(gen) %>% unique() -> genlist
  for (j in 1:nrow(genlist)){
    if (genlist[j,] == 1)next
    for (k in 1:4){
      tmp <- read.delim(paste("/home/owens/working/texanus/linkage_tests/",basename,".",
                            poplist[i],"_",genlist[j,],".",parents[k],".geno.ld.mean",sep=""))
      tmp$location <- poplist[i]
      tmp$gen <- genlist[j,]
      tmp$parent <- parents[k]
      parental_linkage <- rbind(parental_linkage, tmp)
    }
  }
}
parental_linkage %>% filter(location == "BFL",chr=="HanXRQChr17")%>%
  filter(parent == "ann1") %>%
  ggplot(.) + geom_point(aes(x=pos,y=mean_r2)) + facet_wrap(~gen)




