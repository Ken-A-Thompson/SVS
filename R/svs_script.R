## Standing Genetic Variation & Speciation R Analysis script
## Author Ken A. Thompson
### 02 Apr 2017

#load packages
#install.packages('tidyverse')
#install.packages('splitstackshape')
library(cowplot)
library(splitstackshape)
library(tidyverse)

### Run themes

theme_ng1 <- theme(aspect.ratio=1.0,panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   axis.line = element_line(size=1),
                   axis.line.x = element_line(color="black", size = 1),
                   axis.line.y = element_line(color="black", size = 1),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black"),
                   axis.title=element_text(color="black"),
                   axis.title.y=element_text(vjust=0.2, size=12),
                   axis.title.x=element_text(vjust=0.1,size=12),
                   axis.text.x=element_text(size=10),
                   axis.text.y=element_text(size=10),
                   legend.position="none",
                   legend.title=element_blank(),
                   plot.title = element_blank())

theme_fig2 <- theme(aspect.ratio=1.0,panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   axis.line = element_line(size=1),
                   axis.line.x = element_blank(),
                   axis.line.y = element_blank(),
                   axis.ticks=element_blank(),
                   axis.text=element_blank(),
                   axis.title=element_blank(),
                   axis.title.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   legend.position="none",
                   legend.title=element_blank(),
                   plot.title = element_blank())
                   



### End theme


# # Testing script for plotting adaptive walk
# # Generate vectors
# 
# mut.d1 <- c(0,
#             rnorm(20)) %>%
#   abs() %>%
#   sort() %>% 
#   as.vector()
#   
# mut.d2 <- c(0, 
#             rnorm(20)) %>%
#   abs() %>% 
#   sort() %>% 
#   as.vector()
# 
# # Generate a 'hybrid'
# 
# ## Create 'parents'
# anc.d1 <- rep(0, length(mut.d1))
# anc.d2 <- rep(0, length(mut.d2))
# 
# #Generate Hybrids
# muts.d1 <- rev(diff(mut.d1, lag = 1))
# muts.d2 <- rev(diff(mut.d2, lag = 1))
# 
# muts.d1.vec <- c(0, muts.d1)
# muts.d2.vec <- c(0, muts.d2)
# 
# hyb.d1 <- replicate(1000, sum(sample(c(muts.d1.vec, anc.d1), length(muts.d1.vec),replace=F))) #Sum of 'trait value' for D1
# hyb.d2 <- replicate(1000, sum(sample(c(muts.d2.vec, anc.d2), length(muts.d2.vec),replace=F))) #Sum of 'trait value' for D2
# 
# 
## Generate figure for 'adaptive' walk

# 
# df <- data.frame(mut.d1, mut.d2) # Generate data frame for ggplot with both mutation vectors
# 
# ### plot it 
# adapt.walk <- ggplot(df, aes(x = rev(mut.d1), y= rev(mut.d2))) +
#   geom_point() + 
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(rev(mut.d1), n=-1), NA), yend=c(tail(rev(mut.d2), n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open")) + 
#   theme_ng1
# adapt.walk
# 
# #Plot hybrids
# ## Create hybrid dataframe
# 
# hyb.df <- data.frame(hyb.d1, hyb.d2)
# 
# adapt.walk + geom_point(data = hyb.df, aes(x = hyb.d1, y = hyb.d2), alpha = 0.25)


# Import numpy data into R

sgv.phenos <- read.csv("/Users/Ken/Documents/Projects/SVS/data/phenos_K1000_n2_B2_u0.001_alpha0.01.csv", header = F, col.names = 1:ncol(sgv.phenos), check.names = F, na.strings = c("[", "]"))

## Add generation number to dataset
gen <- 1:nrow(sgv.phenos)
sgv.phenos.gen <- cbind(gen, sgv.phenos)

sgv.tidy <- sgv.phenos.gen %>% 
  gather(key = individual, value = pheno, 2:ncol(sgv.phenos.gen)) %>% #convert to long format
  mutate(pheno = gsub("\\[|\\]","", pheno)) %>%  #get rid of square brackets
  mutate(pheno.t = trimws(pheno)) %>% #get rid of leading white space
  mutate(pheno.D = paste(pheno.t, "D", step = "")) %>%  #add 'dummy character' "D" to end of string to enable deletion
  mutate(pheno.t2 = gsub("0.  0. D ","0 0", pheno.D)) %>%  #get rid of double zero oddity
  mutate(pheno.t3 = gsub("D","", pheno.t2)) %>%  #get rid of Dummy
  mutate(pheno.t4 = gsub("  "," ", pheno.t3)) %>%  #get rid of double space
  mutate(pheno.t5 = trimws(pheno.t4)) #get rid of trailing white space

#use cSplit to get X and Y
sgv.tidy.split <- cSplit(sgv.tidy, splitCols = "pheno.t5", sep = " ")

sgv.plot <- sgv.tidy.split %>% 
  select(gen, individual, pheno.t5_1, pheno.t5_2) %>% 
  rename(X = pheno.t5_1) %>% 
  rename(Y = pheno.t5_2) 

write.csv(sgv.plot, file = "R/2017-03-27-tidySVS.csv")

#Finally... now plot it

arrows <- plyr::ddply(na.omit(sgv.plot), ~gen, summarise, meanX = mean(X), meanY= mean(Y))

real.walk.means <- ggplot(arrows[1:20,], aes(x = meanX, y= meanY)) +
  geom_point() + 
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open")) +
  theme_ng1
real.walk.means

real.walk.means + geom_point(data = sgv.plot, aes(x = X, y = Y, colour = gen), alpha = 0.1)

real.walk <- ggplot(na.omit(sgv.plot), aes(x = X, y= Y, colour = gen)) +
 geom_point(aes(colour = gen), alpha = 0.2) + #colour by generation; lighter are later on
 stat_summary(fun.y=mean, geom = "line", aes(group =factor(gen))) +
 labs(x = "Trait 1", y = "Trait 2") +
 theme_ng1
real.walk

# Plot divergent adaptation then hybrids.
## Load in burn-in, both adaptive walks, and hybrid dataframe
sgv.burn <- read.csv()
sgv.adapt1 <- read.csv()
sgv.adapt2 <- read.csv()
sgv.hybrid <- read.csv()

### Fig. 1. SGV with increasing burn-in

sgv.100 <- read.csv('data/ancestor_pop_K1000_n2_B2_u0.001_alpha0.02_gens100_burn.csv', header = F, check.names = F)

sgv.1000 <- read.csv('data/ancestor_pop_K1000_n2_B2_u0.001_alpha0.02_gens1000_burn.csv', header = F, check.names = F)

sgv.5000 <- read.csv('data/ancestor_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_burn.csv', header = F, check.names = F)

sgv.10000 <- read.csv('data/ancestor_pop_K1000_n2_B2_u0.001_alpha0.02_gens10000_burn.csv', header = F, check.names = F)

sgv.20000 <- read.csv('data/ancestor_pop_K1000_n2_B2_u0.001_alpha0.02_gens20000_burn.csv', header = F, check.names = F)

#try distancematrix in hopach
#Within population diversity function returns the mean euclidean distance between individuals in a population
within.pop.diversity <- function(x) {
  dist.mat <- as.matrix(dist(x, method = "euclidean"))
  dist.mat[upper.tri(dist.mat)] <- NA
  pairwise.dist <- as.data.frame(cbind(which(!is.na(dist.mat),arr.ind = TRUE),na.omit(as.vector(dist.mat))))
  pairwise.dist.new <- pairwise.dist %>% 
    mutate(same.ind = row - col) %>% 
    filter(same.ind != 0) %>% 
    select(-same.ind) 
  return(mean(pairwise.dist.new$V3))
}


# Figure 2: Adaptive walks with hybrids from DNM only (A-B) and DNM + SGV

## Dependent functions

## Create function to convert numpy array to long format
numpy.as.long <- function(x) {
  gen <- 1:nrow(x) ## Create variable with 'generation' number
  phenos.gen <- cbind(gen, x) ## Add generation number to numpy array
  phenos.tidy <- phenos.gen %>% 
    gather(key = individual, value = pheno, 2:ncol(sgv.phenos.gen)) %>% #convert to long format
    mutate(pheno = gsub("\\[|\\]","", pheno)) 
  phenos.tidy.split <- cSplit(phenos.tidy, splitCols = "pheno", sep = " ") #use cSplit to get X and Y
  return(phenos.tidy.split)
}

#Create function that joins the two datasets 
join.and.group <- function(x, y){
  data1 <- numpy.as.long(x)
  data2 <- numpy.as.long(y)
  data1$group <- rep("A", nrow(data1))
  data2$group <- rep("B", nrow(data1))
  plot.data <- rbind(data1, data2) %>% 
    rename(X = pheno_1) %>% 
    rename(Y = pheno_2)
  plot.data$group <- as.factor(plot.data$group)
  return(plot.data)
}

## Load data
### DNM only
DNM1pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51_DNM.csv', header = F, col.names = 1:ncol(SGV1pos), check.names = F, na.strings = c("[", "]"))
DNM2pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49_DNM.csv', header = F, col.names = 1:ncol(SGV1pos), check.names = F, na.strings = c("[", "]"))
DNM2neg <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt-0.51--0.51_DNM.csv', header = F, col.names = 1:ncol(SGV1pos), check.names = F, na.strings = c("[", "]"))

###SGV + DNM
SGV1pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51.csv', header = F, col.names = 1:ncol(SGV1pos), check.names = F, na.strings = c("[", "]"))
SGV2pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49.csv', header = F, col.names = 1:ncol(SGV2pos), check.names = F, na.strings = c("[", "]"))
SGV1neg <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt-0.51--0.51.csv', header = F, col.names = 1:ncol(SGV2pos), check.names = F, na.strings = c("[", "]"))


# Fig. 2A Adaptation to parallel optima with only DNM
Parallel.DNM.Data <- join.and.group(DNM1pos, DNM2pos)
Fig2A.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))

Fig.2A.Data <- Parallel.DNM.Data %>% 
  group_by(group, gen) %>% 
  summarise(meanX = mean(X), meanY = mean(Y)) %>% 
  as.data.frame() %>% 
  add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
  add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
  arrange(group, gen)

Fig.2A <- ggplot(Fig.2A.Data, aes(x = meanX, y= meanY, colour = group)) +
  geom_point() + 
  geom_point(data = Fig2A.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2A.Data[1:16,]) +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2A.Data[17:32,]) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  theme_fig2
Fig.2A
 
# Fig. 2B Adaptation to parallel optima with only DNM
Divergent.DNM.Data <- join.and.group(DNM1pos, DNM2neg)
Fig2B.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))

Fig.2B.Data <- Divergent.DNM.Data %>%  #can make this a function a bit later on if need be...
  group_by(group, gen) %>% 
  summarise(meanX = mean(X), meanY = mean(Y)) %>% 
  as.data.frame() %>% 
  add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
  add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
  arrange(group, gen)

Fig.2B <- ggplot(Fig.2B.Data, aes(x = meanX, y= meanY, colour = group)) +
  geom_point() + 
  geom_point(data = Fig2B.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2B.Data[1:16,]) +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2B.Data[17:32,]) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  theme_fig2
Fig.2B

## Fig. 2C - Adaptation to parallel optima from both SGV and DNM
Parallel.SGV.Data <- join.and.group(SGV1pos, SGV2pos)
Fig2C.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))

Fig.2C.Data <- Parallel.SGV.Data %>% 
  group_by(group, gen) %>% 
  summarise(meanX = mean(X), meanY = mean(Y)) %>% 
  as.data.frame() %>% 
  add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
  add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
  arrange(group, gen)
  
Fig.2C <- ggplot(Fig.2C.Data, aes(x = meanX, y= meanY, colour = group)) +
  geom_point() + 
  geom_point(data = Fig2C.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2C.Data[1:16,]) +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2C.Data[17:32,]) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  theme_fig2
Fig.2C

#Fig. 2D Adaptation to divergent optima from both SGV and DNM
Divergent.SGV.Data <- join.and.group(SGV1pos, SGV1neg)
Fig2D.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y"))) # load hybrid data (need to make it so it's unique)

Fig.2D.Data <- Divergent.SGV.Data %>%
  group_by(group, gen) %>% 
  summarise(meanX = mean(X), meanY = mean(Y)) %>% 
  as.data.frame() %>% 
  add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
  add_row(group = "B", gen = 0, meanX = 0, meanY = 0)  %>% #ancestor
  arrange(group, gen)

Fig.2D <- ggplot(Fig.2D.Data, aes(x = meanX, y= meanY, colour = group)) +
  geom_point() + 
  geom_point(data = Fig2D.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
 # geom_point(data = Fig2C.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2D.Data[1:16,]) +
  geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2D.Data[17:32,]) +
  # geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
  #              arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2C.Data[17:32,]) + #plots parallal and divergent this and above)
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  theme_fig2
Fig.2D

#Create multipanel figure
Fig.2 <- plot_grid(Fig.2A, Fig.2B, Fig.2C, Fig.2D, labels = c("A", "B", "C", "D"))




