library(ggplot2)
library(dplyr)
library(reshape)

data <- read.table("/home/owens/working/texanus/texanus.gatk.noQTL.50misscount.5dp.10mac.bi.2.meanQ",colClasses="numeric")
colnames(data) <- c("G1","G2")
labels <- read.delim("/home/owens/working/texanus/texanus.noQTL.samplelist.txt",header=F)
colnames(labels) <- "ID"
info <- read.delim("/home/owens/working/texanus/texanus.sampleinfo.txt",header=T)

data2 <- cbind(data, labels)
data2 <- merge(data2, info)
data2 <- melt(data2, id.var=c("ID","year","gen","plot_type","type","location"))
data2$value <- as.numeric(data2$value)

data2$Name <- factor(data2$Name, 
                     levels=data2[order(data2$Type),]$Name)

data2 %>%
 filter(location == "KPC") %>%
  ggplot(aes(x = ID, y = value, fill = factor(variable)))+
  geom_bar(stat = "identity", width = 1.1)+
  facet_grid(~gen, scales = "free", space = "free",switch="x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        panel.margin = unit(0.1, "lines"), 
        strip.background = element_blank(),
        strip.switch.pad.grid = unit(0.0, "lines"))+
  theme(plot.margin = unit(c(0,0.5,0.5,0.0), "cm")) +
  ylab("q-value")+
  scale_fill_manual(values=c("#123eabff","#ff0000ff"),
                    name="Groups",
                    labels=c("Annuus (1)", "Petiolaris (2)")) +
  theme(strip.text.x = element_text(angle=90,size=6)) -> k2
