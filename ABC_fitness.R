#Function for running an ABC on selection. 
#Input a list of genotype frequencies 
library(abc)
library(tidyverse)
chrlist <- read.delim("/home/owens/working/texanus/rqtl/chrlist.txt",header=F)
popsize <- read_tsv("texanus_popsizes.txt")
abc_out <- tibble(parent=as.character(),fitness=as.numeric(),
                  chr=as.character(),pop=as.character(),
                  start=as.numeric(),type=as.character())

point_estimates <- tibble(parent=as.character(),stat=as.character(),fitness=as.numeric(),
                          chr=as.character(),start=as.numeric(),
                          type=as.character())

seed_bank_percent <- 0.1
chosen_pop <- "BFL"
chosen_start <- 0
chosen_chr <- "HanXRQChr01"
N_sim <- 1000000

for (chosen_chr in chrlist[,1]){
  print(paste("For", chosen_chr))
  for (chosen_pop in poplist){
    print(paste("For",chosen_pop))
    window_frequency %>% filter(chr == chosen_chr,location == chosen_pop) %>% 
      ungroup() %>%dplyr::select(start) %>% unique() %>% .$start -> starts
    
    
    for (chosen_start in starts){
      print(paste("For",chosen_start))
      window_frequency %>% filter(location == chosen_pop, chr == chosen_chr, start == chosen_start) %>% 
        ungroup() %>% dplyr::select(gen) %>% unique() %>% floor() %>% .$gen %>% sort()-> sampled_generations
      #Ginput = (previous generation) ann1, ann2, ann3, deb, generation gap, sample_size
      Ginput <- matrix(NA,nrow=(length(sampled_generations) -1),ncol=6)
      Gfreq_final <- c()
      Ne_input <- matrix(NA, nrow=(length(sampled_generations) -1),ncol=10)
      n_row <- 1
      for (i in 2:length(sampled_generations)){
        window_frequency %>% filter(chr == chosen_chr, 
                                    start == chosen_start, 
                                    location == chosen_pop,
                                    floor(gen) == sampled_generations[i]) %>%
          ungroup() -> window_frequency_tmp
        window_frequency_tmp %>%
          dplyr::select(rel_freq) %>% .$rel_freq -> tmp_current
        
        window_frequency %>% filter(chr == chosen_chr, 
                                    start == chosen_start, 
                                    location == chosen_pop,
                                    floor(gen) == sampled_generations[(i-1)]) %>%
          ungroup() %>%
          dplyr::select(rel_freq) %>% .$rel_freq -> tmp_previous
        for (j in 1:length(tmp_previous)){
          if (tmp_previous[j] == 0){
            tmp_previous[j] = 0.01 #If it's not there, put it there at 1%
          }
        }
        
        window_frequency_tmp  %>%
          head(n=1) %>%
          .$total_alleles -> sample_size
        gen_gap <- sampled_generations[i] - sampled_generations[(i-1)] 
        
        n_spot = 1
        for (j in (sampled_generations[(i-1)]+1):(sampled_generations[(i)])){
          Ne_tmp <- popsize %>% filter(Plot == chosen_pop, Generation == j) %>% .$PopSize
          if (is.na(Ne_tmp)){
            Ne_tmp <- popsize %>% filter(Plot == chosen_pop) %>% 
              summarise(estimated = median(PopSize, na.rm=T)) %>% .$estimated
          }
          Ne_input[n_row,n_spot] <- Ne_tmp *2
          n_spot = n_spot+ 1
        }
        #This is the input for the simulation
        Ginput[n_row,] <- c(tmp_previous, gen_gap,sample_size)
        n_row = n_row+1
        #This is the real measured output, for comparison with simulations.
        Gfreq_final <- c(Gfreq_final, tmp_current)
      }

      oo<-abcloop(N=N_sim,Ne_input=Ne_input,Ginput=Ginput,lb=0,ub=2,sampled_generations=sampled_generations)
      
      aout<-abc(target=Gfreq_final,param=as.matrix(oo[[1]]),sumstat=as.matrix(oo[[2]]),tol=0.0005,method="ridge")

      wout<-matrix(NA,nrow=(N_sim*0.0005),ncol=3)
      wout[,1]<-aout$adj.values[,1] 
      wout[,2]<-aout$adj.values[,2]
      wout[,3]<-aout$adj.values[,3] 
      wout <- as.data.frame(wout)
      colnames(wout) <- c("ann1","ann2","ann3")
      wout <- gather(wout, key) 
      wout %>%
        ggplot(.,aes(y=value,x=1)) +
        geom_violin(aes(colour=key),alpha=0.5) +
        geom_hline(yintercept = 1)
      
      wout$chr <- chosen_chr
      wout$pop <- chosen_pop
      wout$start <- chosen_start
      wout$type <- "full_generations"
      colnames(wout) <- c("parent","fitness","chr","pop","start","type")
      abc_out <- rbind(abc_out, wout)
      
      parents_list <- c("ann1","ann2","ann3");
      values_list <- c("min2.5","median","mean","mode","max97.5")
      value_positions <- c(2:6)
      aout.summary <- summary(aout)
      for (i in 1:3){
        current_parent <- parents_list[i]
        for (j in 1:5){
          point_estimates_tmp <- tibble(parent=current_parent,
                                        stat=values_list[j],
                                        fitness=aout.summary[value_positions[j],i],
                                        chr=chosen_chr,
                                        pop=chosen_pop,
                                        start=chosen_start,
                                        type="full_generations")
          point_estimates <- rbind(point_estimates, point_estimates_tmp)
        }

      }

      
      
      
    }
  }
}

onesim<-function(W,Ginput, Ne_input, generations, sampled_generations){
  
  Gfreq_total <- c()
  ## post selection expected gen freqs
  for (i in 1:nrow(Ginput)){
    if (i > 1){
      Gfreq <- (Ginput[i,1:4] * (1-seed_bank_percent)) + (Ginput[(i-1),1:4] * seed_bank_percent)
    }else{
      Gfreq <- Ginput[i,1:4]
    }
    for (j in 1:Ginput[i,5]){
      Ne <- Ne_input[i,j]
      if (Ne == 0){
        next
      }
      Gfreq <- Gfreq * W
      Gcnts<-rmultinom(n=1,size=Ne,prob=Gfreq)
      Gfreq<-Gcnts/Ne
    }
    N_sampled <- Ginput[i,6]
    Gcnts_sim <- rmultinom(n=1,size=N_sampled,prob=Gfreq)
    Gfreq_sim <- Gcnts_sim/sum(Gcnts_sim)
    Gfreq_total <- c(Gfreq_total, Gfreq_sim)
  }
  return(t(Gfreq_total)) 

}

abcloop<-function(N=1000,Ne_input=Ne_input,Ginput=Ginput,lb=0,ub=2,sampled_generations=sampled_generations){

  simfr<-matrix(NA,nrow=N,ncol=4*(length(sampled_generations) -1)) #This is the result of the ABC, to be checked against real data
  s<-matrix(NA,nrow=N,ncol=4) #This is the parameter space tested.
  for(x in 1:N){
    ## defines the relative fitness of each
    s[x,1:3]<-runif(n=3,min=lb,max=ub)
    s[x,4] <- 1
    W<- s[x,1:4]
    simfr[x,]<-onesim(W=W,Ne_input=Ne_input,Ginput=Ginput,generations=generations,sampled_generations=sampled_generations)
  }    
  out<-list(s,simfr)
  return(out)
}   


pdf("ABC_fitness_window_v3_popsize_10seedbank.pdf")

for (chosen_chr in chrlist[,1]){
  print(
    abc_out %>%
      filter(chr == chosen_chr) %>% 
      ggplot(.,aes(fill=parent,fitness)) +
      geom_density(aes(fill=parent,color=parent),alpha=0.5) +
     # xlim(c(-1,2)) +
      facet_grid(pop~start,scales="free_y") + 
      geom_vline(xintercept=1,color=RColorBrewer::brewer.pal(4,"Set1")[4]) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      ggtitle(paste(chosen_chr,"ABC-fitness estimates. 1 = Debilis fitness")) +
      coord_flip(xlim=c(0,2),ylim=c(0,5))
  )
  print(
    point_estimates %>%
      filter(chr == chosen_chr)  %>%
      group_by(chr, start, pop, parent) %>%
      summarize(mode = fitness[which(stat == "mode")],
                max97.5 = fitness[which(stat == "max97.5")],
                min2.5 = fitness[which(stat == "min2.5")]) %>%
      ggplot(.,aes(x=parent,y=mode,color=parent)) +
      geom_point() +
      geom_errorbar(aes(ymin=min2.5, ymax=max97.5)) +
      geom_hline(yintercept = 1) +
      facet_grid(pop~start) +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      ggtitle(paste(chosen_chr,"ABC-fitness estimates. 1 = Debilis fitness"))
  )
}
dev.off()

point_estimates %>%
  filter(chr == chosen_chr)  %>%
  group_by(chr, start, pop, parent) %>%
  summarize(mode = fitness[which(stat == "mode")],
            max97.5 = fitness[which(stat == "max97.5")],
            min2.5 = fitness[which(stat == "min2.5")]) %>%
  ggplot(.,aes(x=parent,y=mode,color=parent)) +
  geom_point() +
  geom_errorbar(aes(ymin=min2.5, ymax=max97.5)) +
  geom_hline(yintercept = 1) +
  facet_grid(pop~start) +
  scale_color_brewer(palette = "Set1")





deb_values <- tibble(chr=character(),pop=character(),start=numeric(),
                     parent=character(),fitness=numeric(),type=character(),
                     stat=character())
  
for (chosen_chr in chrlist[,1]){
  for (chosen_pop in poplist){
    print(paste("For",chosen_pop))
    window_frequency %>% filter(chr == chosen_chr,location == chosen_pop) %>% 
      ungroup() %>%dplyr::select(start) %>% unique() %>% .$start -> starts
      for (chosen_start in starts){
        tmp <- tibble(chr=as.character(chosen_chr),pop=as.character(chosen_pop),start=as.numeric(chosen_start),
                      parent=as.character("deb"),fitness=as.numeric(1), type=as.character("full_generations"),
                      stat=as.character("mode"))
        deb_values <- rbind(deb_values,tmp)
      }
  }
}
pdf("ABC_fitness_window_popsize_pointestimates_density_10seedbank_v2.pdf",width=6,height=4)
point_estimates %>%
  rbind(., deb_values) %>% 
  filter(stat == "mode") %>%
  group_by(chr, pop, start) %>%
  mutate(max_fitness = max(fitness)) %>%
  mutate(rel_fitness = fitness/max_fitness) %>% 
  ggplot(.) + 
  geom_density(aes(rel_fitness,color=parent,fill=parent),alpha=0.2) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Density") + xlab("Relative fitness") +
  facet_wrap(~pop,nrow=2)
dev.off()
pdf("ABC_fitness_window_popsize_pointestimates_10seedbank_v2.pdf",width=15,height=2.5)
point_estimates %>%
  rbind(., deb_values) %>% 
  filter(stat == "mode") %>%
  group_by(chr, pop, start) %>%
  mutate(max_fitness = max(fitness)) %>%
  mutate(rel_fitness = fitness/max_fitness) %>% 
  mutate(chr_n = str_replace(chr, "HanXRQChr","")) %>% 
  ggplot(.) + 
  geom_segment(aes(x=start,xend=(start+20000000),
                   y=as.factor(parent),yend=as.factor(parent),
                   color=rel_fitness),size=2) +
  facet_grid(pop~chr_n) +
  scale_color_distiller(palette = "RdYlBu", limits=c(0.5,1),name="Relative_fitness") +
  theme_bw() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Parent") +xlab("")
  

point_estimates %>%
  rbind(., deb_values) %>% 
  filter(stat == "mode") %>%
  group_by(chr, pop, start) %>%
  mutate(max_fitness = max(fitness)) %>%
  mutate(rel_fitness = fitness/max_fitness) %>% 
  filter(parent == "deb") %>%
  mutate(chr_n = str_replace(chr, "HanXRQChr","")) %>% 
  ggplot(.) + 
  geom_segment(aes(x=start,xend=(start+20000000),
                   y=as.factor(parent),yend=as.factor(parent),
                   color=rel_fitness),size=2) +
  facet_grid(pop~chr_n) +
  scale_color_distiller(palette = "RdYlBu", limits=c(0.5,1),name="Relative_fitness") +
  theme_bw() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Parent") +xlab("")
dev.off()



densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  return(td$x[maxDens])
}

####ARCHIVES
abc_out_noseedbank <- abc_out
abc_out_10seedbank <- abc_out
  
abc_out_noseedbank %>% group_by(chr, pop, start, parent) %>% 
  summarise(fitness = densMode(fitness)) %>% 
  ungroup() %>%
  rbind(., deb_values) %>% 
  group_by(chr, pop, start) %>%
  mutate(max_fitness = max(fitness)) %>%
  mutate(rel_fitness_noseed = fitness/max_fitness)  -> abc_out_noseedbank_pointestimate

abc_out_10seedbank %>% group_by(chr, pop, start, parent) %>% 
  summarise(fitness = densMode(fitness)) %>% 
  ungroup() %>%
  rbind(., deb_values) %>% 
  group_by(chr, pop, start) %>%
  mutate(max_fitness = max(fitness)) %>%
  mutate(rel_fitness_10seed = fitness/max_fitness) -> abc_out_10seedbank_pointestimate

pdf("Comparing_seedbank_fitness.v0.pdf")
merge(abc_out_noseedbank_pointestimate %>% select(-fitness, -max_fitness), 
      abc_out_10seedbank_pointestimate %>% select(-fitness, -max_fitness)) %>% 
  ggplot(.,aes(x=rel_fitness_10seed, y=rel_fitness_noseed)) +
  geom_point() + geom_smooth(method=lm) + facet_wrap(~parent,nrow=2)
dev.off()

#Plotting population sizes
pdf("Texanus_population_sizes.pdf",height=4,width=4)
popsize %>% 
  filter(!is.na(PopSize)) %>%
  ggplot(.) + 
  geom_line(aes(x=Generation,y=PopSize),na.rm=T) + 
  facet_wrap(~Plot,scale="free") + 
  theme_bw() +
  ylab("Population size")
dev.off()





