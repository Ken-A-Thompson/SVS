## Standing Genetic Variation & Speciation R Analysis script
## Author Ken A. Thompson
### 02 Apr 2017

#load packages
# install.packages('cowplot')
# install.packages('flexclust')
# install.packages('pbapply')
# install.packages('tidyverse')
# install.packages('splitstackshape')
library(cowplot)
library(data.table)
library(flexclust)
# library(pbapply)
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
                   # legend.position="none",
                   # legend.title=element_blank(),
                   plot.title = element_blank())
                   
### End themes

### Helper functions

# Create function that calculates pairwise distance
between.pop.diversity <- function(x, y){
  ##first thing's first: delete 'locus 1' which is not useful information here because it is an automatic fixed difference that is a neutral mutation.
  x %>% select(-1)
  y %>% select(-1)
  #Block of code to append loci with all zeroes to 'sister' population
  ## first population
  zeroes.1 <- matrix(0L, nrow(x), ncol(y))
  homologs.1 <- as.data.frame(zeroes.1)
  pop.w.homologs.1 <- cbind(x, homologs.1)
  ## second population
  zeroes.2 <- matrix(0L, nrow(y), ncol(x))
  homologs.2 <- as.data.frame(zeroes.2)
  pop.w.homologs.2 <- cbind(homologs.2, y) #do this in opposite order as above so that loci are 'aligned'
  ## Get the pairwise euclidean distances between the matrices
  dist.mat <- as.matrix(dist2(pop.w.homologs.1, pop.w.homologs.2, method = "euclidean")) #dist2 calculates for all pairs of rows across datasets
  dist.mat[upper.tri(dist.mat)] <- NA
  pairwise.dist <- as.data.frame(cbind(which(!is.na(dist.mat),arr.ind = TRUE),na.omit(as.vector(dist.mat))))
  return(mean(pairwise.dist$V3))
}

#Mean hybrid fitness for two traits at optimum of 0.5
hybrid.fitness.n2.opt0.5 <- function(x,y){ 
  euclid <- (sqrt((x - 0.5)^2 + (y - 0.5)^2))^2
  fitness <- exp(-1 * euclid)
  return(mean(fitness))
}

#Create function that joins the two 'adaptive walk' datasets 
join.and.group <- function(x, y){
  data1 <- numpy.as.long(x)
  data2 <- numpy.as.long(y)
  data1$group <- rep("A", nrow(data1)) #this and below assigns each a unique group identifier
  data2$group <- rep("B", nrow(data1))
  plot.data <- rbind(data1, data2) %>% #combine the two datasets
    rename(X = pheno_1) %>%  #rename the phenotypes to x and Y
    rename(Y = pheno_2)
  plot.data$group <- as.factor(plot.data$group) # make sure 'group' is coded as a factor
  return(plot.data)
}

mutation.matrices <- function(x) { #this function works with lapply to take all the list items (generations) and make them into 'nice' matrices. 
  cleaned.matrix <- t(as.data.frame(x))
  fin.matrix <- data.frame(cleaned.matrix) %>% 
    gather(key = individual, value = locus, 1:ncol(cleaned.matrix)) %>% #convert to long format
    mutate(locus = gsub("\\[|\\]","", locus)) %>% #lose the square brackets, they suck!
    cSplit(splitCols = "locus", sep = " ") %>%  #use cSplit to get all loci in unique columns
    select(-individual) 
  return(fin.matrix)
}


# Function to convert numpy array to long format
numpy.as.long <- function(x) {
  gen <- 1:nrow(x) ## Create variable with 'generation' number
  phenos.gen <- cbind(gen, x) ## Add generation number to numpy array
  phenos.tidy <- phenos.gen %>% 
    gather(key = individual, value = pheno, 2:ncol(x)) %>% #convert to long format
    mutate(pheno = gsub("\\[|\\]","", pheno)) 
  phenos.tidy.split <- cSplit(phenos.tidy, splitCols = "pheno", sep = " ") #use cSplit to get X and Y
  return(phenos.tidy.split)
}

# Function to remove columns that have only zeros in them (modified from SO to work in list loops)
remove_zero_cols <- function(x) {
  df <- as.data.frame(x)
  rem_vec <- NULL
  for(i in 1:ncol(df)){
    this_sum <- summary(df[,i])
    zero_test <- length(which(this_sum == 0))
    if(zero_test == 6) {
      rem_vec[i] <- names(df)[i]
    }
  }
  features_to_remove <- rem_vec[!is.na(rem_vec)]
  rem_ind <- which(names(df) %in% features_to_remove)
  df <- df[,-rem_ind]
  return(df)
}

#Within population diversity function returns the mean euclidean distance between individuals (rows) in a population
within.pop.diversity <- function(x) {
  dist.mat <- as.matrix(dist(sample_n(x, 1000), method = "euclidean")) #have the sample_n in there for now to limit to 1000 individuals (makes it run faster)
  dist.mat[upper.tri(dist.mat)] <- NA
  pairwise.dist <- as.data.frame(cbind(which(!is.na(dist.mat),arr.ind = TRUE),na.omit(as.vector(dist.mat))))
  pairwise.dist.new <- pairwise.dist %>% 
    mutate(same.ind = row - col) %>% 
    filter(same.ind != 0) %>% 
    select(-same.ind) 
  return(mean(pairwise.dist.new$V3))
}

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

# Hybrid + parent fitness calculation
## Using equation of Fraisse et al. 2016



hybrid.fitness.n2.opt0.5(Fig2C.Hybrids$X, Fig2C.Hybrids$Y)

hybrid.fitness.n2.opt0.5(Fig2A.Hybrids$X, Fig2A.Hybrids$Y)

###############################################
################### FIGURES ###################
###############################################

# Fig. 1. SGV with increasing burn-in
## Load base data for fig (burn-in population sampled at different time points)

Fig1.RawData <- fread('data/ancestor_pop_K10000_n2_B2_u0.001_alpha0.02_gens2000_burn.csv', header = F, check.names = F, na.strings = c("[", "]")) #use fread bc read.csv breaks R

# Make each generation a list item
Fig1.List <- as.list(data.frame(t(Fig1.RawData)))

# Use the 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance (takes about a minute)
Fig1.CleanList <- lapply(Fig1.List, mutation.matrices) 

# Calculate within population diversity for every generation (can take a while, but currently am sampling just 1000 individuals from the 10000-long data)
Fig1.SGV <- data.frame(sapply(Fig1.CleanList, within.pop.diversity))

#Make data for Fig 1A
Fig1A.Data <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(Fig1.SGV))), c("Generation", "SGV")) #Create dataframe
Fig1A.Data$SGV <- Fig1.SGV$sapply.Fig1.CleanList..within.pop.diversity. #Bring in data from sapply above with all within-pop pairwise distances
Fig1A.Data$Generation <- as.vector(1:nrow(Fig1.SGV) * 100) #The second number is the frequency with which the data were sampled
Fig.1A.Data.Sorted <- Fig1A.Data[order(Fig1B.Data$Generation),]

#Create dataset to point at 'burn-in' arrow
Fig1A.Arrow <- data.frame(x1 = last(Fig.1A.Data.Sorted$Generation), x2 = last(Fig.1A.Data.Sorted$Generation), y1 = 0.65*last(Fig.1A.Data.Sorted$SGV), y2 = 0.8*last(Fig.1A.Data.Sorted$SGV))

#Plot figure 1A
Fig.1A <- ggplot(Fig.1A.Data.Sorted, aes(x = Generation, y = SGV)) +
  geom_point() + 
  labs(x = "Generation #",
       y = "Within-population genetic diversity") +
  geom_smooth() + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), arrow=arrow(length=unit(0.3,"cm"), type = "open"), size = 1, data = Fig1A.Arrow) +
  theme_ng1
Fig.1A

#Plot Figure 1B
## number of mutations in the population

# Need to count length of mutations (non-empty) in each generation... plot that
# Remove columns that have only zeros; these mutations were present before selection but not after (i.e., they were just lost)
Muts.No.Zeros <- lapply(Fig1.CleanList, remove_zero_cols)

#Count ncols of each 'generation'
Muts.ColCnt <- t(data.frame(lapply(Muts.No.Zeros, ncol)))

#Make data for Fig. 2B
#Make data for Fig 1A
Fig1B.Data <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(Fig1.SGV))), c("Generation", "Nmuts")) #Create dataframe
Fig1B.Data$Nmuts <- Muts.ColCnt[,1] #Bring in data from sapply above with all within-pop pairwise distances
Fig1B.Data$Generation <- as.vector(1:nrow(Fig1.SGV) * 100) #The second number is the frequency with which the data were sampled
Fig1B.Data <- add_row(Fig1B.Data, Generation = 0, Nmuts = 1) #Root at 1 (start with one)
Fig.1B.Data.Sorted <- Fig1B.Data[order(Fig1B.Data$Generation),]

#Create dataset to point at 'burn-in' arrow
Fig1B.Arrow <- data.frame(x1 = last(Fig.1B.Data.Sorted$Generation), x2 = last(Fig.1B.Data.Sorted$Generation), y1 = 0.65*last(Fig.1B.Data.Sorted$Nmuts), y2 = 0.8*last(Fig.1B.Data.Sorted$Nmuts))

Fig.1B <- ggplot(Fig.1B.Data.Sorted, aes(x = Generation, y = Nmuts)) +
  geom_point() + 
  labs(x = "Generation #",
       y = "Number of mutations in population") +
  geom_smooth() + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), arrow=arrow(length=unit(0.3,"cm"), type = "open"), size = 1, data = Fig1B.Arrow) +
  theme_ng1
Fig.1B


## Create multi-panel figure for Fig 1

Fig.1 <- plot_grid(Fig.1A, Fig.1B, labels = c("a", "b"))


# Figure 2: Adaptive walks with hybrids from DNM only (A-B) and DNM + SGV

## Load data
### DNM only
DNM1pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51_DNM.csv', header = F, check.names = F, na.strings = c("[", "]"))
DNM2pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49_DNM.csv', header = F, check.names = F, na.strings = c("[", "]"))
DNM2neg <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt-0.51--0.51_DNM.csv', header = F, check.names = F, na.strings = c("[", "]"))

###SGV + DNM
SGV1pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
SGV2pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49.csv', header = F, check.names = F, na.strings = c("[", "]"))
SGV1neg <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt-0.51--0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))

## Euclidean distance between populations
### Load in the mutation matrices
SGV1pos.muts <- read.csv('data/parent_pop_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
SGV2pos.muts <- read.csv('data/parent_pop_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49.csv', header = F, check.names = F, na.strings = c("[", "]"))



# Make each generation an item in a list
matrix.list.SGV1pos.muts <- as.list(data.frame(t(SGV1pos.muts)))
matrix.list.SGV2pos.muts <- as.list(data.frame(t(SGV2pos.muts)))

#use the 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
clean.list.SGV1pos.muts <- lapply(matrix.list.SGV1pos.muts, mutation.matrices) #next get everything into a 'good' list.
clean.list.SGV2pos.muts <- lapply(matrix.list.SGV2pos.muts, mutation.matrices) #next get everything into a 'good' list.

# Now that we have a list item with a clean matrix representing each population's 
# mutations for each array, need to get a column that contains their pairwise euclidean
# distance at each generation. I.e., 'dist2' for each pair of list items in corresponding dataset



Fig.2A.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.SGV1pos.muts, y =  clean.list.SGV2pos.muts)) #use 'mapply' to iterate the calculation over between all elements of the lists. 

# also need to fix the hybid code because it always loads the most recent one

# Fig. 2A Adaptation to parallel optima with only DNM
Parallel.DNM.Data <- join.and.group(DNM1pos, DNM2pos)
Fig2A.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))

# 'Genomic divergence'
Fig2A.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.SGV1pos.muts, y =  clean.list.SGV2pos.muts)) #use 'mapply' to iterate the calculation over between all elements of the lists. 

## Create data for Fig 2A
###Should make these functions if there's time
Fig.2A.Data <- Parallel.DNM.Data %>%
  group_by(group, gen) %>% 
  summarise(meanX = mean(X), meanY = mean(Y)) %>% 
  as.data.frame() %>% 
  mutate(divergence = c(Fig2A.Divergence[,1], Fig2A.Divergence[,1])) %>% 
  add_row(group = "A", gen = 0, meanX = 0, meanY = 0, divergence = 0) %>% #ancestor
  add_row(group = "B", gen = 0, meanX = 0, meanY = 0, divergence = 0) %>%  #ancestor
  arrange(group, gen)

Fig.2A <- ggplot(Fig.2A.Data, aes(x = meanX, y= meanY, colour = divergence)) +
  # scale_colour_gradientn(colours = rainbow(7)) +
  scale_colour_gradient(low = 'red', high = 'blue') +
  geom_point() + 
  geom_point(data = Fig2A.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.1) +
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




