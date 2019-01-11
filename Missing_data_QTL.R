#Plotting missing data in mapping population

library(tidyverse)
parents <-c("ann1","ann2","ann3","deb")

cross_data <- tibble(locus=character(),chr=character(),cm=numeric(),genotype=numeric(),parent=character())
for (parent in parents){
tmp <- read_csv(paste("/home/owens/working/texanus/texanus.gatk.GBS1.BC1.maf10.bi.GBS2",parent,"sites.recode.csvr",sep=""),col_names=F)
tmp <- tmp[26:nrow(tmp),]
tmp <- tmp %>% gather(colnames(tmp)[4:length(colnames(tmp))], key = "sample", value = "genotype") 
colnames(tmp) <- c("locus","chr","cm","sample","genotype")
tmp <- tmp %>% filter(genotype != "-")
tmp$parent <- parent
cross_data <- rbind(cross_data,tmp)
}

pdf("Missing_data_all.pdf",width=10,height=10)
all_chr <- unique(sort(cross_data$chr))
for (chosen_parent in parents){
  for (i in all_chr){
    print(
      cross_data%>% 
        filter(parent == chosen_parent, chr == i) %>%
        mutate(factor_cm = as.factor(cm), factor_sample = as.factor(sample)) %>%
        group_by(sample) %>% mutate(ancestry = sum(genotype == "A")/sum(!is.na(genotype))) %>% 
        ggplot(.,aes(x=factor_cm,y=fct_reorder(factor_sample, ancestry),fill=genotype)) +
        geom_tile() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        ggtitle(paste(chosen_parent," BC1, missing data, Chr ",i,sep=" ")) +
        ylab("Samples") + xlab("Loci")
    )
  }
}
dev.off()


cross_data_imputed <- tibble(locus=character(),chr=character(),cm=numeric(),genotype=numeric(),parent=character())
for (parent in parents){
  tmp <- read_csv(paste("/home/owens/working/texanus/texanus.gatk.GBS1.BC1.maf10.bi.GBS2",parent,"sites.recode.imputed.csvr",sep=""),col_names=F)
  tmp <- tmp[26:nrow(tmp),]
  tmp <- tmp %>% gather(colnames(tmp)[4:length(colnames(tmp))], key = "sample", value = "genotype") 
  colnames(tmp) <- c("locus","chr","cm","sample","genotype")
  tmp <- tmp %>% filter(genotype != "-")
  tmp$parent <- parent
  cross_data_imputed <- rbind(cross_data_imputed,tmp)
}

cross_data$type <- "real"
cross_data_imputed$type <- "imputed"
cross_data_full <- rbind(cross_data,cross_data_imputed)

pdf("Imputed_data_all.pdf",width=10,height=10)
all_chr <- unique(sort(cross_data$chr))
for (chosen_parent in parents){
  for (i in all_chr){
    print(
      cross_data_full%>% 
        filter(parent == chosen_parent, chr == i) %>%
        mutate(factor_cm = as.factor(cm), factor_sample = as.factor(sample)) %>%
        group_by(sample) %>% mutate(ancestry = sum(genotype == "A")/sum(!is.na(genotype))) %>% 
        filter(type == "imputed") %>% 
        ggplot(.,aes(x=factor_cm,y=fct_reorder(factor_sample, ancestry))) +
        geom_tile(fill="grey") +
        geom_tile(data=cross_data_full %>%
                    filter(parent == chosen_parent, chr == i) %>%
                    mutate(factor_cm = as.factor(cm), factor_sample = as.factor(sample)) %>%
                    group_by(sample) %>% mutate(ancestry = sum(genotype == "A")/sum(!is.na(genotype))) %>%
                    filter(type == "real"), aes(fill=genotype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        ggtitle(paste(chosen_parent," BC1, missing data, Chr ",i,sep=" ")) +
        ylab("Samples") + xlab("Loci")
    )
  }
}
dev.off()

cross_data_combined <- tibble(locus=character(),chr=character(),cm=numeric(),genotype=numeric(),parent=character())
for (parent in parents[c(2,4)]){
  tmp <- read_csv(paste("/home/owens/working/texanus/texanus.gatk.GBS1.BC1.maf10.bi.GBS2",parent,"sites.recode.combined.csvr",sep=""),col_names=F)
  tmp <- tmp[26:nrow(tmp),]
  tmp <- tmp %>% gather(colnames(tmp)[4:length(colnames(tmp))], key = "sample", value = "genotype") 
  colnames(tmp) <- c("locus","chr","cm","sample","genotype")
  tmp <- tmp %>% filter(genotype != "-")
  tmp$parent <- parent
  cross_data_combined <- rbind(cross_data_combined,tmp)
}

pdf("combined_data_all.pdf",width=10,height=10)
all_chr <- unique(sort(cross_data$chr))
for (chosen_parent in parents[c(2,4)]){
  for (i in all_chr){
    print(
      cross_data_combined%>% 
        filter(parent == chosen_parent, chr == i) %>%
        mutate(factor_cm = as.factor(cm), factor_sample = as.factor(sample)) %>%
        group_by(sample) %>% mutate(ancestry = sum(genotype == "A")/sum(!is.na(genotype))) %>% 
        ggplot(.,aes(x=factor_cm,y=fct_reorder(factor_sample, ancestry),fill=genotype)) +
        geom_tile() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        ggtitle(paste(chosen_parent," BC1, missing data, Chr ",i,sep=" ")) +
        ylab("Samples") + xlab("cM")
    )
  }
}
dev.off()


