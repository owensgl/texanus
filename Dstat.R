library(tidyverse)
library(ggmap)
library(viridis)
library(gridExtra)
directory <- "/home/owens/working/texanus/WGS"
info <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv")
pop_loc <- read_tsv("pop_loc_allnum.txt") %>% rename(population = pop)
usa <- map_data('state')
states <- map_data("state")
target_state <- subset(states, region %in% c("texas"))
lat_range <- c(25, 34)
long_range <- c(-106,-96)



deb_alleles <- read_tsv(paste(directory,"/Texanus.tranche90.snp.gwas.90.bi.Deballeles.90.5.txt",sep=""),
                        col_names = c("name","deb_alleles"))

inner_join(deb_alleles,info) %>% full_join(.,pop_loc) %>% 
  filter(taxon == "annuus") %>% mutate(type = case_when(lat < 33 & long > -106 ~ "Tex",
                                                        TRUE ~"Cen")) %>%
  filter(type == "Cen") %>% summarize(max_cen = max(deb_alleles), min_cen = min(deb_alleles),
                                      mean_cen = mean(deb_alleles)) -> cen_range

plot1 <- inner_join(deb_alleles,info) %>% full_join(.,pop_loc) %>% 
  filter(taxon == "annuus" | taxon== "argophyllus") %>% 
  mutate(type = case_when(lat < 33 & long > -106 ~ "Tex",
                          TRUE ~"Cen")) %>%  
  filter(type == "Tex") %>% 
  mutate(population = case_when(taxon == "argophyllus" ~ "ARG",
                                TRUE ~ population)) %>%
  mutate(long = case_when(taxon == "argophyllus" ~ -95,
                          TRUE ~ long)) %>%
  ggplot(.,aes(x=fct_reorder(population,long), y=deb_alleles)) + geom_boxplot(aes()) +
  theme_bw() +
  geom_hline(yintercept=cen_range$max_cen,color="blue") +
  geom_hline(yintercept=cen_range$min_cen,color="blue") +
  geom_hline(yintercept=cen_range$mean_cen,linetype="dotted",color="blue") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("Ann.Tex population by longitude") +
  ylab("Proportion of debilis alleles present") 





plot2 <- ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_point(data=inner_join(deb_alleles,info) %>% 
                    full_join(.,pop_loc) %>% 
                    filter(taxon == "annuus") %>% 
                    mutate(type = case_when(lat < 33 & long > -106 ~ "Tex",
                                            TRUE ~"Cen")) %>%  
                    filter(type == "Tex") %>%
               group_by(population,lat,long) %>%
               summarize(mean_deb = mean(deb_alleles)),
                  aes(x=long, y=lat,color=mean_deb,size=3)) +
  scale_color_viridis(name="Mean proportion of\ndebilis alleles\npresent") +theme_bw() +
  xlab("Longitude") + ylab("Latitude") + 
  coord_map(xlim = long_range,ylim = lat_range)


####Now looking at Dstat

Deb_Dstat <- read_tsv(paste(directory,"/Texanus.tranche90.snp.gwas.ANN-DEB.Dstat.txt",sep=""))

plot3 <- Deb_Dstat %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  ggplot(aes(x=fct_reorder(name,long),y=Dstat)) + geom_violin() +
  theme_bw() +
  geom_hline(yintercept = 0,color="red") +
  ylab("D") + xlab("Ann.Tex samples by longitude") +
  ggtitle("(Ann.Tex Ann.Cen Deb) ABBA-BABA") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plot4 <- Deb_Dstat %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  group_by(name,population,lat,long) %>%
  summarise(mean_Dstat = mean(Dstat,na.rm=T)) %>% 
  filter(!is.na(mean_Dstat)) %>%
  ggplot(aes(x=fct_reorder(population,long),y=mean_Dstat)) + geom_violin() +
  theme_bw() +
  geom_hline(yintercept = 0,color="red") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab("D") + xlab("Ann.Tex population by longitude") +
  ggtitle("(Ann.Tex Ann.Cen Deb) ABBA-BABA") 


Arg_Dstat <- read_tsv(paste(directory,"/Texanus.tranche90.snp.gwas.ANN-ARG.Dstat.txt",sep=""))

plot5 <- Arg_Dstat %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  ggplot(aes(x=fct_reorder(name,long),y=Dstat)) + geom_violin() +
  theme_bw() +
  geom_hline(yintercept = 0,color="red") +
  ylab("D") + xlab("Ann.Tex samples by longitude") +
  ggtitle("(Ann.Tex Ann.Cen Arg) ABBA-BABA") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

plot6 <- Arg_Dstat %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  group_by(name,population,lat,long) %>%
  summarise(mean_Dstat = mean(Dstat,na.rm=T)) %>% 
  filter(!is.na(mean_Dstat)) %>%
  ggplot(aes(x=fct_reorder(population,long),y=mean_Dstat)) + geom_violin() +
  theme_bw() +
  geom_hline(yintercept = 0,color="red") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab("D") + xlab("Ann.Tex population by longitude") +
  ggtitle("(Ann.Tex Ann.Cen Arg) ABBA-BABA") 

plot7 <- Arg_Dstat %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  group_by(name,population,lat,long) %>%
  summarise(mean_Dstat = mean(Dstat,na.rm=T)) %>% 
  filter(!is.na(mean_Dstat)) %>% inner_join(.,deb_alleles) %>% 
  ggplot(.,aes(x=mean_Dstat,y=deb_alleles)) + geom_point(aes(color=population)) +
  theme_bw() +
  geom_hline(yintercept=cen_range$max_cen) +
  geom_hline(yintercept=cen_range$min_cen) +
  geom_hline(yintercept=cen_range$mean_cen,linetype="dotted") +
  ylab("Proportion of debilis alleles present") + 
  xlab("ARG D") +
  ggtitle("The relationship between deb alleles\nand Arg introgression")


pdf("Texanus.tranche90.snp.gwas.90.bi.Deballeles.abba.v1.pdf",height=16,width=8)
grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7, nrow = 4)
dev.off()

plot8 <- Deb_Dstat %>%
  rename(name = P2) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  filter(long < -116) %>% 
  rename(P2 = name) %>%
  select(P1, P2, P3, Dstat) %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  ggplot(aes(x=fct_reorder(name,long),y=Dstat)) + geom_violin() +
  theme_bw() +
  geom_hline(yintercept = 0,color="red") +
  ylab("D") + xlab("Ann.Tex samples by longitude") +
  ggtitle("(Ann.Tex Ann.Cal Deb) ABBA-BABA") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plot9 <-Deb_Dstat %>%
  rename(name = P2) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  filter(long < -116) %>% 
  rename(P2 = name) %>%
  select(P1, P2, P3, Dstat) %>%
  rename(name = P1) %>%
  inner_join(.,info) %>% full_join(.,pop_loc) %>% 
  group_by(name,population,lat,long) %>%
  summarise(mean_Dstat = mean(Dstat,na.rm=T)) %>% 
  filter(!is.na(mean_Dstat)) %>%
  ggplot(aes(x=fct_reorder(population,long),y=mean_Dstat)) + geom_violin() +
  theme_bw() +
  geom_hline(yintercept = 0,color="red") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab("D") + xlab("Ann.Tex population by longitude") +
  ggtitle("(Ann.Tex Ann.Cal Deb) ABBA-BABA") 

pdf("Texanus.tranche90.snp.gwas.90.bi.Cal.abba.v1.pdf",height=8,width=8)
grid.arrange(plot3, plot4, plot8,plot9, nrow = 2)
dev.off()

