## Standing Genetic Variation & Speciation R Analysis script
## Author Ken A. Thompson
### 2018-01-08

#load packages
# install.packages('cowplot')
# install.packages('flexclust')
# install.packages('pbapply')
# install.packages('tidyverse')
# install.packages('splitstackshape')
# library(cowplot)
# library(data.table)
# library(flexclust)
# # library(pbapply)
# library(splitstackshape)
library(akima)
library(fields)
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
                   # legend.position="none",
                   # legend.title=element_blank(),
                   plot.title = element_blank())
# 
# theme_fig2 <- theme(aspect.ratio=1.0,panel.background = element_blank(),
#                    panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    panel.border=element_blank(),
#                    axis.line = element_line(size=1),
#                    axis.line.x = element_blank(),
#                    axis.line.y = element_blank(),
#                    axis.ticks=element_blank(),
#                    axis.text=element_blank(),
#                    axis.title=element_blank(),
#                    axis.title.y=element_blank(),
#                    axis.title.x=element_blank(),
#                    axis.text.x=element_blank(),
#                    axis.text.y=element_blank(),
#                    # legend.position="none",
#                    # legend.title=element_blank(),
#                    plot.title = element_blank())
# 
# theme_Fig4 <- theme(aspect.ratio=1.0,panel.background = element_blank(),
#                    panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    panel.border=element_blank(),
#                    axis.line = element_line(size=1),
#                    axis.line.x = element_line(color="black", size = 1),
#                    axis.line.y = element_line(color="black", size = 1),
#                    axis.ticks=element_line(color="black"),
#                    axis.text=element_text(color="black"),
#                    axis.title=element_text(color="black"),
#                    axis.title.y=element_text(vjust=0.2, size=12),
#                    axis.title.x=element_text(vjust=0.1,size=12),
#                    axis.text.x=element_text(size=10),
#                    axis.text.y=element_text(size=10),
#                    legend.position="none",
#                    legend.title=element_blank(),
#                    plot.title = element_blank())
#                    
# End themes
################
### Figure 2 ###
################

### Hybrid load vs. number of mutations for both Divergent and parallel selection

# read data
Fig2.D1 <- read.csv('SVS_plots/2018-01-07-SVS_PlotsWCSVs/data2.csv')
Fig2.D2 <- read.csv('SVS_plots/2018-01-07-SVS_PlotsWCSVs/data3.csv')

#join (temp until I create them together and optimize this)
Fig2.Data <- rbind(Fig2.D1, Fig2.D2)

# seperate into three datasets for optima
Fig2.Data.Far <-
  Fig2.Data %>% 
  filter(opt == 0.75) %>% 
  group_by(nmuts, regime) %>% 
  mutate(hl.mean = mean(hybrid.load), hl.sd = sd(hybrid.load))

Fig.2C <- ggplot(Fig2.Data.Far, aes(x = nmuts, y = hl.mean, colour = regime)) +
  geom_point() + 
  labs(x = "Number of ancestral mutations",
       y = "Hybrid load") +
  geom_smooth(method = "loess") +
  geom_errorbar(aes(ymin = hl.mean - hl.sd, ymax = hl.mean + hl.sd)) + 
  ylim(0, 0.085) +
  theme_ng1
Fig.2C

################
### Figure 3 ###
################

### Hybrid load across environments for different amts of mutation

Fig3.Data <- read.csv('SVS_plots/2018-01-08-SVS_Plots.txt/data3.csv')

Fig.Data.PlotReady <-
  Fig3.Data %>% 
  group_by(angle, nmuts) %>% 
  mutate(hl.mean = mean(hybrid.load), hl.sd = sd(hybrid.load))

Fig.3 <- ggplot(Fig.Data.PlotReady, aes(x = angle, y = hl.mean, colour = factor(nmuts))) +
  geom_point() + 
  labs(x = "Angle separating environments",
       y = "Hybrid load") +
  geom_smooth(method = "loess") +
  # geom_errorbar(aes(ymin = hl.mean - hl.sd, ymax = hl.mean + hl.sd)) + 
  # ylim(0, 0.085) +
  # geom_vline(xintercept = 90) + 
  theme_ng1
Fig.3

################
### Figure 4 ###
################

### Hybrid load effects on fitness across environments

Fig4.Data <- read_csv('data/d1_heatmap_percentile.csv', col_names = c("nmuts", "hybrid.load", "angle", "meanfit", "maxfit", "ninetypercent"))


Fig4A.Data <- interp(Fig4.Data$angle,Fig4.Data$hybrid.load,Fig4.Data$meanfit)
Fig4A.Data <-  image.plot(Fig4A.Data)
abline(v = 90)

Fig4B.Data <- interp(Fig4.Data$angle,Fig4.Data$hybrid.load,Fig4.Data$ninetypercent)
Fig4B <- image.plot(Fig4B2.Data)

pdf('fig4a.pdf',width=4,height=3) 


??fields
# ## Create mult
# 
# 
# ### Plot it
# 
# Fig.3F <- ggplot(Fig.3F.Data, aes(x = 100*Gen, y= Dist, colour = Group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Genetic divergence") +
#   geom_smooth() + 
#   xlim(0,5000) +
#   ylim(0, 38) +
#   theme_ng1
# Fig.3F
# 
# ## Create mult
# 
# 
# 
# 
# ### Helper functions
# 
# # Create function that calculates pairwise distance
# between.pop.diversity <- function(x, y){
#   ##first thing's first: delete 'locus 1' which is not useful information here because it is an automatic fixed difference that is a neutral mutation.
#   x %>% select(-1)
#   y %>% select(-1)
#   #Block of code to append loci with all zeroes to 'sister' population
#   ## first population
#   zeroes.1 <- matrix(0L, nrow(x), ncol(y))
#   homologs.1 <- as.data.frame(zeroes.1)
#   pop.w.homologs.1 <- cbind(x, homologs.1)
#   ## second population
#   zeroes.2 <- matrix(0L, nrow(y), ncol(x))
#   homologs.2 <- as.data.frame(zeroes.2)
#   pop.w.homologs.2 <- cbind(homologs.2, y) #do this in opposite order as above so that loci are 'aligned'
#   ## Get the pairwise euclidean distances between the matrices
#   dist.mat <- as.matrix(dist2(pop.w.homologs.1, pop.w.homologs.2, method = "manhattan")) #dist2 calculates for all pairs of rows across datasets
#   dist.mat[upper.tri(dist.mat)] <- NA
#   pairwise.dist <- as.data.frame(cbind(which(!is.na(dist.mat),arr.ind = TRUE),na.omit(as.vector(dist.mat))))
#   return(mean(pairwise.dist$V3))
# }
# 
# # Function that calculates pairwise distance between populatins that share SGV; need ancestral columns to be defined as 'sgv.end'
# ## Now works in the rare case that a pop does not have any mutations beyond the ancestral ones
# between.pop.diversity.SGV <- function(x, y, sgv.end){
#   ##define a range of columns that includes all ancestral SGV
#   x.only.anc <- x %>% select(1:sgv.end) 
#   y.only.anc <- y  %>% select(1:sgv.end)
#   ## define the range of columns that are de novo mutation by dropping all sgv
#   x.no.anc <- x  %>% select(-(1:sgv.end))
#   y.no.anc <- y  %>% select(-(1:sgv.end))
#   #Block of code to append loci with all zeroes to 'sister' population
#   ## first population
#   zeroes.1 <- matrix(0L, nrow(x.only.anc), ncol(y.no.anc))
#   homologs.1 <- as.data.frame(zeroes.1)
#   pop.w.homologs.1 <- x.only.anc
#   if(ncol(x.no.anc) > 0){
#     pop.w.homologs.1 <- cbind(pop.w.homologs.1, x.no.anc)
#   }
#   if(ncol(homologs.1) > 0){
#     pop.w.homologs.1 <- cbind(pop.w.homologs.1, homologs.1)
#   }
#   ## second population
#   zeroes.2 <- matrix(0L, nrow(y.only.anc), ncol(x.no.anc))
#   homologs.2 <- as.data.frame(zeroes.2)
#   pop.w.homologs.2 <- y.only.anc
#   if(ncol(homologs.2) > 0){
#     pop.w.homologs.2 <- cbind(pop.w.homologs.2, homologs.2)
#   }
#   if(ncol(y.no.anc) > 0){
#     pop.w.homologs.2 <- cbind(pop.w.homologs.2, y.no.anc)
#   }
#   ## Get the pairwise euclidean distances between the matrices
#   dist.mat <- as.matrix(dist2(pop.w.homologs.1, pop.w.homologs.2, method = "manhattan")) #dist2 calculates for all pairs of rows across datasets
#   dist.mat[upper.tri(dist.mat)] <- NA
#   pairwise.dist <- as.data.frame(cbind(which(!is.na(dist.mat),arr.ind = TRUE),na.omit(as.vector(dist.mat))))
#   return(mean(pairwise.dist$V3))
# }
# 
# 
# #Mean hybrid fitness for two traits at optimum of 0.5
# hybrid.fitness.n2.opt <- function(x,y, opt, sigma.seln = 1){ 
#   euclid <- (sqrt((x - opt)^2 + (y - opt)^2))^2
#   fitness <- exp(-1 * sigma.seln * euclid)
#   return(fitness)
# }
# 
# #Create function that joins the two 'adaptive walk' datasets 
# join.and.group <- function(x, y){
#   data1 <- numpy.as.long(x)
#   data2 <- numpy.as.long(y)
#   data1$group <- rep("A", nrow(data1)) #this and below assigns each a unique group identifier
#   data2$group <- rep("B", nrow(data1))
#   plot.data <- rbind(data1, data2) %>% #combine the two datasets
#     rename(X = pheno_1) %>%  #rename the phenotypes to x and Y
#     rename(Y = pheno_2)
#   plot.data$group <- as.factor(plot.data$group) # make sure 'group' is coded as a factor
#   return(plot.data)
# }
# 
# mutation.matrices <- function(x) { #this function works with lapply to take all the list items (generations) and make them into 'nice' matrices. 
#   cleaned.matrix <- t(as.data.frame(x))
#   fin.matrix <- data.frame(cleaned.matrix) %>% 
#     gather(key = individual, value = locus, 1:ncol(cleaned.matrix)) %>% #convert to long format
#     mutate(locus = gsub("\\[|\\]","", locus)) %>% #lose the square brackets, they suck!
#     cSplit(splitCols = "locus", sep = " ") %>%  #use cSplit to get all loci in unique columns
#     select(-individual) 
#   return(fin.matrix)
# }
# 
# # Function to convert numpy array to long format
# numpy.as.long <- function(x) {
#   gen <- 1:nrow(x) ## Create variable with 'generation' number
#   phenos.gen <- cbind(gen, x) ## Add generation number to numpy array
#   phenos.tidy <- phenos.gen %>% 
#     gather(key = individual, value = pheno, 2:ncol(x)) %>% #convert to long format
#     mutate(pheno = gsub("\\[|\\]","", pheno)) 
#   phenos.tidy.split <- cSplit(phenos.tidy, splitCols = "pheno", sep = " ") #use cSplit to get X and Y
#   return(phenos.tidy.split)
# }
# 
# # Function to summarize large population array into mean phenotype by generation
# phenotype.summary <- function(x) {
#   phenotype.data <- x %>%
#     group_by(group, gen) %>% 
#     summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#     as.data.frame() %>% 
#     add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% # starting phenotype
#     add_row(group = "B", gen = 0, meanX = 0, meanY = 0)  %>% # starting phenotype
#     arrange(group, gen) %>% 
#     rowwise() %>% 
#     mutate(TraitMean = mean(c(meanX, meanY), na.rm=T))  #calculate a 'trait mean' of X and Y.
#   return(phenotype.data)
# }
# 
# # Function to remove columns that have only zeros in them (modified from SO to work in list loops)
# remove_zero_cols <- function(x) {
#   df <- as.data.frame(x)
#   rem_vec <- NULL
#   for(i in 1:ncol(df)){
#     this_sum <- summary(df[,i])
#     zero_test <- length(which(this_sum == 0))
#     if(zero_test == 6) {
#       rem_vec[i] <- names(df)[i]
#     }
#   }
#   features_to_remove <- rem_vec[!is.na(rem_vec)]
#   rem_ind <- which(names(df) %in% features_to_remove)
#   df <- df[,-rem_ind]
#   return(df)
# }
# 
# #can probably delete this block when Matt's revised code works.
# 
# # need to get rid of the 'df' column if there is one... figure it out!
# 
# # Modify remove_zero_cols for preserving ancestral SGV but removing DNM
# remove_zero_cols.sgv <- function(x) {
#   dnm.cols <- x[,dnm.start:ncol(x)]
#   sgv.cols <- x[,1:sgv.end]
#   df <- as.data.frame(dnm.cols)
#   rem_vec <- NULL
#   for(i in 1:ncol(df)){
#     this_sum <- summary(df[,i])
#     zero_test <- length(which(this_sum == 0))
#     if(zero_test == 6) {
#       rem_vec[i] <- names(df)[i]
#     }
#   }
#   features_to_remove <- rem_vec[!is.na(rem_vec)]
#   rem_ind <- which(names(df) %in% features_to_remove)
#   df <- df[,-rem_ind]
#   joint.data <- cbind(sgv.cols, df)
#   if('df' %in% colnames(joint.data)) {
#     joint.data <- select(joint.data, -df)
#   }
#   return(joint.data)
# }
# 
# #if there is a column called 'df' then just give me sgvcols
# #if there is no column called df then give me the cbind.
# 
# 
# 
# 
# #Within population diversity function returns the mean euclidean distance between individuals (rows) in a population
# within.pop.diversity <- function(x) {
#   dist.mat <- as.matrix(dist(sample_n(x, 1000), method = "manhattan")) #have the sample_n in there for now to limit to 1000 individuals (makes it run faster)
#   dist.mat[upper.tri(dist.mat)] <- NA
#   pairwise.dist <- as.data.frame(cbind(which(!is.na(dist.mat),arr.ind = TRUE),na.omit(as.vector(dist.mat))))
#   pairwise.dist.new <- pairwise.dist %>% 
#     mutate(same.ind = row - col) %>% 
#     filter(same.ind != 0) %>% 
#     select(-same.ind) 
#   return(mean(pairwise.dist.new$V3))
# }
# 
# ###############################################
# ################### FIGURES ###################
# ###############################################
# 
# # Fig. 1. SGV with increasing burn-in
# ## Load base data for fig (burn-in population sampled at different time points)
# 
# Fig1.RawData <- fread('data/ancestor_pop_K10000_n2_B2_u0.001_alpha0.02_gens5000_burn_rep1_sgv.csv', header = F, check.names = F, na.strings = c("[", "]")) #use fread bc read.csv breaks R
# 
# # Make each generation a list item
# Fig1.List <- as.list(data.frame(t(Fig1.RawData)))
# 
# # Use the 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance (takes about a minute)
# Fig1.CleanList <- lapply(Fig1.List, mutation.matrices) 
# 
# # Calculate within population diversity for every generation (can take a while, but currently am sampling just 1000 individuals from the 10000-long data)
# Fig1.SGV <- data.frame(sapply(Fig1.CleanList, within.pop.diversity))
# 
# #Make data for Fig 1A
# Fig1A.Data <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(Fig1.SGV))), c("Generation", "SGV")) #Create dataframe
# Fig1A.Data$SGV <- Fig1.SGV$sapply.Fig1.CleanList..within.pop.diversity. #Bring in data from sapply above with all within-pop pairwise distances
# Fig1A.Data$Generation <- as.vector(1:nrow(Fig1.SGV) * 100) #The second number is the frequency with which the data were sampled
# Fig.1A.Data.Sorted <- Fig1A.Data[order(Fig1A.Data$Generation),]
# 
# #Create datasets to point at 'burn-in' arrow
# Fig1A.Vert.Arrow <- data.frame(x1 = last(Fig.1A.Data.Sorted$Generation), x2 = last(Fig.1A.Data.Sorted$Generation), y1 = 0.65*last(Fig.1A.Data.Sorted$SGV), y2 = 0.8*last(Fig.1A.Data.Sorted$SGV))
# 
# Fig1A.Hor.Arrow <- data.frame(x1 = 0.35*last(Fig.1A.Data.Sorted$Generation), x2 = 0.2*last(Fig.1A.Data.Sorted$Generation), y1 = first(Fig.1A.Data.Sorted$SGV), y2 = first(Fig.1A.Data.Sorted$SGV))
# 
# ## Print ancestor less often if need be
# 
# #Plot figure 1A
# Fig.1A <- ggplot(Fig.1A.Data.Sorted, aes(x = Generation, y = SGV)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Within-population genetic diversity") +
#   geom_smooth() + 
#   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), arrow=arrow(length=unit(0.3,"cm"), type = "open"), size = 1, data = Fig1A.Vert.Arrow) +
#   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), arrow=arrow(length=unit(0.3,"cm"), type = "open"), size = 1, data = Fig1A.Hor.Arrow) +
#   theme_ng1
# Fig.1A
# 
# #Plot Figure 1B
# ## number of mutations in the population
# 
# # Need to count length of mutations (non-empty) in each generation... plot that
# # Remove columns that have only zeros; these mutations were present before selection but not after (i.e., they were just lost)
# # Choose more founders; higher carrying capacity.
# 
# Muts.No.Zeros <- lapply(Fig1.CleanList, remove_zero_cols)
# 
# #Count ncols of each 'generation'
# Muts.ColCnt <- t(data.frame(lapply(Muts.No.Zeros, ncol)))
# 
# #Make data for Fig. 1B
# Fig1B.Data <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(Fig1.SGV))), c("Generation", "Nmuts")) #Create dataframe
# Fig1B.Data$Nmuts <- Muts.ColCnt[,1] #Bring in data from sapply above with all within-pop pairwise distances
# Fig1B.Data$Generation <- as.vector(1:nrow(Fig1.SGV) * 100) #The second number is the frequency with which the data were sampled
# Fig1B.Data <- add_row(Fig1B.Data, Generation = 0, Nmuts = 1) #Root at 1 (start with one)
# Fig.1B.Data.Sorted <- Fig1B.Data[order(Fig1B.Data$Generation),]
# 
# #Create dataset to point at 'burn-in' arrow
# Fig1B.Vert.Arrow <- data.frame(x1 = last(Fig.1B.Data.Sorted$Generation), x2 = last(Fig.1B.Data.Sorted$Generation), y1 = 0.65*last(Fig.1B.Data.Sorted$Nmuts), y2 = 0.8*last(Fig.1B.Data.Sorted$Nmuts))
# Fig1B.Hor.Arrow <- data.frame(x1 = 0.35*last(Fig.1B.Data.Sorted$Generation), x2 = 0.2*last(Fig.1B.Data.Sorted$Generation), y1 = first(Fig.1B.Data.Sorted$Nmuts), y2 = first(Fig.1B.Data.Sorted$Nmuts))
# 
# Fig.1B <- ggplot(Fig.1B.Data.Sorted, aes(x = Generation, y = Nmuts)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Number of mutations in population") +
#   geom_smooth() + 
#   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), arrow=arrow(length=unit(0.3,"cm"), type = "open"), size = 1, data = Fig1B.Vert.Arrow) +
#   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), arrow=arrow(length=unit(0.3,"cm"), type = "open"), size = 1, data = Fig1A.Hor.Arrow) +
#   theme_ng1
# Fig.1B
# 
# 
# ## Create multi-panel figure for Fig 1
# 
# Fig.1 <- plot_grid(Fig.1A, Fig.1B, labels = c("a", "b"))
# 
# #############
# ### FIG 2 ###
# #############
# 
# # Figure 2: Adaptive walks with hybrids from DNM only (A-B) and DNM + SGV
# 
# ## Load in phenotypes at every sampled generation
# ### DNM only
# DNM1pos <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# DNM2pos <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.249-0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# DNM2neg <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.249--0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# ###SGV + DNM
# SGV1pos <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.251-0.251_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# SGV2pos <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.249-0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# SGV1neg <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.249--0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# 
# # ## Load in the mutation matrices
# # ### DNM only
# # 
# # DNM1pos.muts <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # DNM2pos.muts <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.249-0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # DNM2neg.muts <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt-0.249--0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # 
# # ###SGV + DNM
# # SGV1pos.muts <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.251-0.251_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # SGV2pos.muts <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.249-0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # SGV2neg.muts <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt-0.249--0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # 
# # ## Make each generation an item in a list
# # matrix.list.DNM1pos.muts <- as.list(data.frame(t(DNM1pos.muts)))
# # matrix.list.DNM2pos.muts <- as.list(data.frame(t(DNM2pos.muts)))
# # matrix.list.DNM2neg.muts <- as.list(data.frame(t(DNM2neg.muts)))
# # 
# # matrix.list.SGV1pos.muts <- as.list(data.frame(t(SGV1pos.muts)))
# # matrix.list.SGV2pos.muts <- as.list(data.frame(t(SGV2pos.muts)))
# # matrix.list.SGV2neg.muts <- as.list(data.frame(t(SGV2neg.muts)))
# # 
# # ## Use the 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# # ## Warnings are OK
# # clean.list.DNM1pos.muts <- lapply(matrix.list.DNM1pos.muts, mutation.matrices)
# # clean.list.DNM2pos.muts <- lapply(matrix.list.DNM2pos.muts, mutation.matrices)
# # clean.list.DNM2neg.muts <- lapply(matrix.list.DNM2neg.muts, mutation.matrices)
# # 
# # ## Use the 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# # ## Warnings are OK
# # clean.list.SGV1pos.muts <- lapply(matrix.list.SGV1pos.muts, mutation.matrices)
# # clean.list.SGV2pos.muts <- lapply(matrix.list.SGV2pos.muts, mutation.matrices)
# # clean.list.SGV2neg.muts <- lapply(matrix.list.SGV2neg.muts, mutation.matrices)
# # 
# # ### Now that we have a list item with a clean matrix representing each population's 
# # ### mutations for each array, get a column that contains their pairwise euclidean
# # ### distance at each generation for each pair of list items in corresponding dataset
# # 
# # ### DNM, parallel
# # Fig.2A.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.DNM1pos.muts, y =  clean.list.DNM2pos.muts))
# # ### DNM, divergent
# # Fig.2B.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.DNM1pos.muts, y =  clean.list.DNM2neg.muts))
# # ### SGV, parallel
# # Fig.2C.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.SGV1pos.muts, y =  clean.list.SGV2pos.muts))
# # ### SGV, divergent
# # Fig.2D.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.SGV1pos.muts, y =  clean.list.SGV2neg.muts))
# 
# #Load all hybrids
# Fig2A.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_0.249-0.249_opt2_0.251-0.251_dnm.csv',header = F, col.names = c("X", "Y")))
# 
# Fig2C.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_0.249-0.249_opt2_0.251-0.251_sgv.csv',header = F, col.names = c("X", "Y")))
# 
# 
# # Fig. 2A Adaptation to parallel optima with only DNM
# Parallel.DNM.Data <- join.and.group(DNM1pos, DNM2pos)
# 
# ## Create data for Fig 2A
# Fig.2A.Data.Full <- Parallel.DNM.Data %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   # mutate(divergence = c(Fig.2A.Divergence[,1], Fig.2A.Divergence[,1])) %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
#   arrange(group, gen)
# 
# Fig.2A.Data.1 <- Fig.2A.Data.Full[1:51,]
# Fig.2A.Data.2 <- Fig.2A.Data.Full[52:102,]
# 
# Fig.2A.Data.Trim.1 <- Fig.2A.Data.1[seq(1, nrow(Fig.2A.Data.1), 3),] 
# Fig.2A.Data.Trim.2 <- Fig.2A.Data.2[seq(1, nrow(Fig.2A.Data.2), 3),] 
# 
# Fig.2A.Data <- rbind(Fig.2A.Data.Trim.1, Fig.2A.Data.Trim.2)
# 
# #Define fitness contours using well-adapted populations
# 
# Neg.Contour <- expand.grid(X = -1 * Fig2A.Hybrids$X, Y = -1 * Fig2C.Hybrids$Y)
# Neg.Contour$Z = hybrid.fitness.n2.opt(Neg.Contour$X, Neg.Contour$Y, opt = -0.25, sigma.seln = 20)
# 
# Pos.Contour <- expand.grid(X = Fig2A.Hybrids$X, Y = Fig2A.Hybrids$Y)
# Pos.Contour$Z = hybrid.fitness.n2.opt(Pos.Contour$X, Pos.Contour$Y, opt = 0.25, sigma.seln = 20)
# 
# ## Fig 2A
# Fig.2A <- ggplot(Fig.2A.Data, aes(x = meanX, y= meanY, colour = group)) +
#   #scale_colour_gradient2(low = 'red', mid = 'purple', midpoint = 20, high = 'blue', data = Fig.2A.Data[1:17,]) +
#   geom_point() + 
#   geom_point(data = Fig2A.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.1) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   # geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#   #              arrow=arrow(length=unit(0.3,"cm"), type = "closed"), colour = "black", size = 0.7, data = Fig.2A.Data[1:17,]) +
#   # geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#   #              arrow=arrow(length=unit(0.3,"cm"), type = "closed"), colour = "black", size = 0.7, data = Fig.2A.Data[18:34,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), size = 0.4, data = Fig.2A.Data[1:17,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), size = 0.4, data = Fig.2A.Data[18:34,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   geom_contour(data = Pos.Contour, aes(x = X, y = Y, z = Z, colour = Z), breaks=c(0.95, 0.85, 0.75), alpha = 0.5) +
#   theme_fig2
# Fig.2A
# 
# Fig2A.Poster <- Fig.2A + theme(legend.position = "none")
# 
# # Fig. 2B Adaptation to divergent optima with only DNM
# #########
# 
# Divergent.DNM.Data <- join.and.group(DNM1pos, DNM2neg)
# Fig2B.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_-0.249--0.249_opt2_0.251-0.251_dnm.csv',header = F, col.names = c("X", "Y")))
# 
# Fig.2B.Data.Full <- Divergent.DNM.Data %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
# #  mutate(divergence = c(Fig.2B.Divergence[,1], Fig.2B.Divergence[,1])) %>% 
# 
# # Figure 2: Adaptive walks with hybrids from DNM only (A-B) and DNM + SGV
# 
# ## Load data
# ### DNM only
# DNM1pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51_DNM.csv', header = F, check.names = F, na.strings = c("[", "]"))
# DNM2pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49_DNM.csv', header = F, check.names = F, na.strings = c("[", "]"))
# DNM2neg <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt-0.51--0.51_DNM.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# ###SGV + DNM
# SGV1pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
# SGV2pos <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49.csv', header = F, check.names = F, na.strings = c("[", "]"))
# SGV1neg <- read.csv('data/parent_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt-0.51--0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# ## Euclidean distance between populations
# ### Load in the mutation matrices
# SGV1pos.muts <- read.csv('data/parent_pop_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.51-0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
# SGV2pos.muts <- read.csv('data/parent_pop_K10000_n2_B2_u0.001_alpha0.02_gens1500_opt0.49-0.49.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# 
# 
# # Make each generation an item in a list
# matrix.list.SGV1pos.muts <- as.list(data.frame(t(SGV1pos.muts)))
# matrix.list.SGV2pos.muts <- as.list(data.frame(t(SGV2pos.muts)))
# 
# #use the 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# clean.list.SGV1pos.muts <- lapply(matrix.list.SGV1pos.muts, mutation.matrices) #next get everything into a 'good' list.
# clean.list.SGV2pos.muts <- lapply(matrix.list.SGV2pos.muts, mutation.matrices) #next get everything into a 'good' list.
# 
# # Now that we have a list item with a clean matrix representing each population's 
# # mutations for each array, need to get a column that contains their pairwise euclidean
# # distance at each generation. I.e., 'dist2' for each pair of list items in corresponding dataset
# 
# Fig.2A.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.SGV1pos.muts, y =  clean.list.SGV2pos.muts)) #use 'mapply' to iterate the calculation over between all elements of the lists. 
# 
# # also need to fix the hybid code because it always loads the most recent one
# 
# # Fig. 2A Adaptation to parallel optima with only DNM
# Parallel.DNM.Data <- join.and.group(DNM1pos, DNM2pos)
# Fig2A.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))
# 
# # 'Genomic divergence'
# Fig2A.Divergence <- data.frame(mapply(between.pop.diversity, x = clean.list.SGV1pos.muts, y =  clean.list.SGV2pos.muts)) #use 'mapply' to iterate the calculation over between all elements of the lists. 
# 
# ## Create data for Fig 2A
# ###Should make these functions if there's time
# Fig.2A.Data <- Parallel.DNM.Data %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   mutate(divergence = c(Fig2A.Divergence[,1], Fig2A.Divergence[,1])) %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0, divergence = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0, divergence = 0) %>%  #ancestor
#   arrange(group, gen)
# 
# Fig.2A <- ggplot(Fig.2A.Data, aes(x = meanX, y= meanY, colour = divergence)) +
#   # scale_colour_gradientn(colours = rainbow(7)) +
#   scale_colour_gradient(low = 'red', high = 'blue') +
#   geom_point() + 
#   geom_point(data = Fig2A.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.1) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2A.Data[1:16,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2A.Data[17:32,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   theme_fig2
# Fig.2A
#  
# # Fig. 2B Adaptation to parallel optima with only DNM
# Divergent.DNM.Data <- join.and.group(DNM1pos, DNM2neg)
# Fig2B.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))
# 
# Fig.2B.Data <- Divergent.DNM.Data %>%  #can make this a function a bit later on if need be...
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
#   arrange(group, gen)
# 
# Fig.2B.Data.1 <- Fig.2B.Data.Full[1:51,]
# Fig.2B.Data.2 <- Fig.2B.Data.Full[52:102,]
# 
# Fig.2B.Data.Trim.1 <- Fig.2B.Data.1[seq(1, nrow(Fig.2B.Data.1), 3),] 
# Fig.2B.Data.Trim.2 <- Fig.2B.Data.2[seq(1, nrow(Fig.2B.Data.2), 3),] 
# 
# Fig.2B.Data <- rbind(Fig.2B.Data.Trim.1, Fig.2B.Data.Trim.2)
# 
# 
# Fig.2B <- ggplot(Fig.2B.Data, aes(x = meanX, y= meanY, colour = group)) +
#   #scale_colour_gradient2(low = 'red', midpoint = 20, mid = 'purple', high = 'blue') +
#   geom_point() + 
#   geom_point(data = Fig2B.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.1) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), data = Fig.2B.Data[1:17,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), data = Fig.2B.Data[18:34,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   geom_contour(data = Pos.Contour, aes(x = X, y = Y, z = Z, colour = Z), breaks=c(0.95, 0.85, 0.75), alpha = 0.5) +
#   geom_contour(data = Neg.Contour, aes(x = X, y = Y, z = Z, colour = Z), breaks=c(0.95, 0.85, 0.75), alpha = 0.5) +
#   theme_fig2
# Fig.2B
# 
# Fig2B.Poster <- Fig.2B + theme(legend.position = "none")
# 
# ## Fig. 2C - Adaptation to parallel optima from both SGV and DNM
# Parallel.SGV.Data <- join.and.group(SGV1pos, SGV2pos)
# 
# Fig.2C.Data.Full <- Parallel.SGV.Data %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
# #  mutate(divergence = c(Fig.2C.Divergence[,1], Fig.2C.Divergence[,1])) %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
#   arrange(group, gen)
# 
# Fig.2C.Data.1 <- Fig.2C.Data.Full[1:51,]
# Fig.2C.Data.2 <- Fig.2C.Data.Full[52:102,]
# 
# Fig.2C.Data.Trim.1 <- Fig.2C.Data.1[seq(1, nrow(Fig.2C.Data.1), 3),] 
# Fig.2C.Data.Trim.2 <- Fig.2C.Data.2[seq(1, nrow(Fig.2C.Data.2), 3),] 
# 
# Fig.2C.Data <- rbind(Fig.2C.Data.Trim.1, Fig.2C.Data.Trim.2)
# 
# ## Fig 2C
# Fig.2C <- ggplot(Fig.2C.Data, aes(x = meanX, y= meanY, colour = group)) +
#   #scale_colour_gradient2(low = 'red', midpoint = 10, mid = 'purple', high = 'blue') +
#   geom_point() + 
#   geom_point(data = Fig2C.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.1) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), data = Fig.2C.Data[1:17,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), data = Fig.2C.Data[18:34,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   geom_contour(data = Pos.Contour, aes(x = X, y = Y, z = Z, colour = Z), breaks=c(0.95, 0.85, 0.75), alpha = 0.5) +
#   theme_fig2
# Fig.2C
# 
# Fig2C.Poster <- Fig.2C + theme(legend.position = "none")
# 
# #Fig. 2D Adaptation to divergent optima from both SGV and DNM
# Divergent.SGV.Data <- join.and.group(SGV1pos, SGV1neg)
# Fig2D.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_-0.249--0.249_opt2_0.251-0.251_sgv.csv',header = F, col.names = c("X", "Y"))) # load hybrid data
# 
# Fig.2D.Data.Full <- Divergent.SGV.Data %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
# #  mutate(divergence = c(Fig.2D.Divergence[,1], Fig.2D.Divergence[,1])) %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
#   arrange(group, gen)
# 
# Fig.2D.Data.1 <- Fig.2D.Data.Full[1:51,]
# Fig.2D.Data.2 <- Fig.2D.Data.Full[52:102,]
# 
# Fig.2D.Data.Trim.1 <- Fig.2D.Data.1[seq(1, nrow(Fig.2D.Data.1), 3),] 
# Fig.2D.Data.Trim.2 <- Fig.2D.Data.2[seq(1, nrow(Fig.2D.Data.2), 3),] 
# 
# Fig.2D.Data <- rbind(Fig.2D.Data.Trim.1, Fig.2D.Data.Trim.2)
# 
# Fig.2D <- ggplot(Fig.2D.Data, aes(x = meanX, y= meanY, colour = group)) +
#   #scale_colour_gradient2(low = 'red', mid = "magenta", midpoint = 20, high = 'dodgerblue4') +
#   geom_point() + 
#   geom_point(data = Fig2D.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.1) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed") , data = Fig.2D.Data[1:17,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "closed"), data = Fig.2D.Data[18:34,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   geom_contour(data = Pos.Contour, aes(x = X, y = Y, z = Z, colour = Z), breaks=c(0.95, 0.85, 0.75), alpha = 0.5) +
#   geom_contour(data = Neg.Contour, aes(x = X, y = Y, z = Z, colour = Z), breaks=c(0.95, 0.85, 0.75), alpha = 0.5) +
#   theme_fig2
# Fig.2D
# 
# Fig2D.Poster <- Fig.2D + theme(legend.position = "none")
# 
# #Create multipanel figure
# 
# Fig.2 <- plot_grid(Fig2A.Poster, Fig2B.Poster, Fig2C.Poster, Fig2D.Poster, labels = c("A", "B", "C", "D"))
# 
# #Save figures for poster
# ggsave(filename = '/Users/Ken/Dropbox/!Ph.D./!SVS_Standing-Variation-Speciation/SVS_Poster/SVS_poster_figs/Fig1A.pdf' , plot = Fig2A.Poster, scale = )
# ggsave(filename = '/Users/Ken/Dropbox/!Ph.D./!SVS_Standing-Variation-Speciation/SVS_Poster/SVS_poster_figs/Fig1B.pdf' , plot = Fig2B.Poster)
# ggsave(filename = '/Users/Ken/Dropbox/!Ph.D./!SVS_Standing-Variation-Speciation/SVS_Poster/SVS_poster_figs/Fig1C.pdf' , plot = Fig2C.Poster)
# ggsave(filename = '/Users/Ken/Dropbox/!Ph.D./!SVS_Standing-Variation-Speciation/SVS_Poster/SVS_poster_figs/Fig1D.pdf' , plot = Fig2D.Poster)
# Fig.2B <- ggplot(Fig.2B.Data, aes(x = meanX, y= meanY, colour = group)) +
#   geom_point() + 
#   geom_point(data = Fig2B.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2B.Data[1:16,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2B.Data[17:32,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   theme_fig2
# Fig.2B
# 
# ## Fig. 2C - Adaptation to parallel optima from both SGV and DNM
# Parallel.SGV.Data <- join.and.group(SGV1pos, SGV2pos)
# Fig2C.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y")))
# 
# Fig.2C.Data <- Parallel.SGV.Data %>% 
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0) %>%  #ancestor
#   arrange(group, gen)
#   
# Fig.2C <- ggplot(Fig.2C.Data, aes(x = meanX, y= meanY, colour = group)) +
#   geom_point() + 
#   geom_point(data = Fig2C.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2C.Data[1:16,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2C.Data[17:32,]) +
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   theme_fig2
# Fig.2C
# 
# #Fig. 2D Adaptation to divergent optima from both SGV and DNM
# Divergent.SGV.Data <- join.and.group(SGV1pos, SGV1neg)
# Fig2D.Hybrids <- as.data.frame(read.csv('data/hybrid_phenos_K10000_n2_B2_u0.001_alpha0.02_gens1500.csv',header = F, col.names = c("X", "Y"))) # load hybrid data (need to make it so it's unique)
# 
# Fig.2D.Data <- Divergent.SGV.Data %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% #ancestor
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0)  %>% #ancestor
#   arrange(group, gen)
# 
# Fig.2D <- ggplot(Fig.2D.Data, aes(x = meanX, y= meanY, colour = group)) +
#   geom_point() + 
#   geom_point(data = Fig2D.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
#  # geom_point(data = Fig2C.Hybrids, aes(x = X, y = Y, colour = NULL), alpha = 0.2) +
#   labs(x = "Trait 1", y = "Trait 2") +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2D.Data[1:16,]) +
#   geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#                arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2D.Data[17:32,]) +
#   # geom_segment(aes(xend=c(tail(meanX, n=-1), NA), yend=c(tail(meanY, n=-1), NA)),
#   #              arrow=arrow(length=unit(0.3,"cm"), type = "open"), data = Fig.2C.Data[17:32,]) + #plots parallal and divergent this and above)
#   xlim(-0.6,0.6) +
#   ylim(-0.6,0.6) +
#   theme_fig2
# Fig.2D
# 
# #Create multipanel figure
# Fig.2 <- plot_grid(Fig.2A, Fig.2B, Fig.2C, Fig.2D, labels = c("A", "B", "C", "D"))
# 
# #################################
# ##Figure 3: Genomic parallelism##
# #################################
# 
# #########
# ##FIG3A##
# #########
# 
# ### Phenotypes
# Fig3.Pos1.SGVDNM <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.251-0.251_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.SGVDNM <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.249-0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos1.SGVDNM <- read.csv('data/parent_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt0.251-0.251_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.SGVDNM <- read.csv('data/parent_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt0.249-0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Use 'join and group' function to put them into the same database.
# Fig3.SGVDNM.Parallel.Phenos <- join.and.group(Fig3.Pos1.SGVDNM, Fig3.Pos2.SGVDNM)
# 
# #Create summary dataset of each group's phenotypic evolution
# Fig.3A.Data <- Fig3.SGVDNM.Parallel.Phenos %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% # starting phenotype
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0)  %>% # starting phenotype
#   arrange(group, gen) %>% 
#   rowwise() %>% 
#   mutate(TraitMean = mean(c(meanX, meanY), na.rm=T))  #calculate a 'trait mean' of X and Y.
# 
# #Plot each 'adaptive walk' on X and Y as a smooth regression, versus GEN.
# #Fig.3A Parallel
# Fig.3A <- ggplot(Fig.3A.Data, aes(x = 100*gen, y= TraitMean, colour = group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Phenotype") +
#   geom_smooth() + 
#   geom_hline(aes(yintercept = 0.25), colour = "black", linetype = "dashed") + 
#   theme_ng1
# Fig.3A
# 
# #########
# ##FIG3B##
# #########
# 
# #Fig.3B Divergent
# #Load in divergent data
# Fig3.Neg1.SGVDNM <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.249--0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Neg1.SGVDNM <- read.csv('data/parent_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.249--0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# #Join the two datasets being compared
# Fig3B.SGVDNM.Divergent.Phenos <- join.and.group(Fig3.Pos1.SGVDNM, Fig3.Neg1.SGVDNM)
# 
# Fig.3B.Data <- Fig3B.SGVDNM.Divergent.Phenos %>%
#   group_by(group, gen) %>% 
#   summarise(meanX = mean(X), meanY = mean(Y)) %>% 
#   as.data.frame() %>% 
#   add_row(group = "A", gen = 0, meanX = 0, meanY = 0) %>% # starting phenotype
#   add_row(group = "B", gen = 0, meanX = 0, meanY = 0)  %>% # starting phenotype
#   arrange(group, gen) %>% 
#   rowwise() %>% 
#   mutate(TraitMean = mean(c(meanX, meanY), na.rm=T))  #calculate a 'trait mean' of X and Y.
# 
# Fig.3B <- ggplot(Fig.3B.Data, aes(x = 100*gen, y= TraitMean, colour = group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Phenotype") +
#   geom_smooth() + 
#   geom_hline(aes(yintercept = 0.25), colour = "black", linetype = "dashed") + 
#   geom_hline(aes(yintercept = -0.25), colour = "black", linetype = "dashed") + 
#   theme_ng1
# Fig.3B
# 
# #########
# ##FIG3C##
# #########
# 
# ### 
# Fig3.Pos1.DNM <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.DNM <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt0.249-0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Use 'join and group' function to put them into the same database.
# Fig3C.SGV.Parallel.Phenos <- join.and.group(Fig3.Pos1.DNM, Fig3.Pos2.DNM)
# Fig3C.Pos1.DNM <- read.csv('data/parent_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3C.Pos2.DNM <- read.csv('data/parent_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt0.249-0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Use 'join and group' function to put them into the same database.
# Fig3C.SGV.Parallel.Phenos <- join.and.group(Fig3C.Pos1.DNM, Fig3C.Pos2.DNM)
# 
# #Create summary dataset of each group's phenotypic evolution
# Fig.3C.Data <- phenotype.summary(Fig3C.SGV.Parallel.Phenos)
# 
# #Plot each 'adaptive walk' on X and Y as a smooth regression, versus GEN.
# #Fig.3A Parallel
# Fig.3C <- ggplot(Fig.3C.Data, aes(x = 100*gen, y= TraitMean, colour = group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Phenotype") +
#   geom_smooth() + 
#   geom_hline(aes(yintercept = 0.25), colour = "black", linetype = "dashed") + 
#   theme_ng1
# Fig.3C
# 
# #########
# ##FIG3D##
# #########
# 
# ### 
# Fig3.Neg1.DNM <- read.csv('data/parent_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.249--0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Use 'join and group' function to put them into the same database.
# Fig3D.DNM.Divergent.Phenos <- join.and.group(Fig3.Pos1.DNM, Fig3.Neg1.DNM)
# Fig3D.Neg1.DNM <- read.csv('data/parent_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.249--0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Use 'join and group' function to put them into the same database.
# Fig3D.DNM.Divergent.Phenos <- join.and.group(Fig3D.Pos1.DNM, Fig3D.Neg1.DNM)
# 
# #Create summary dataset of each group's phenotypic evolution
# Fig.3D.Data <- phenotype.summary(Fig3D.DNM.Divergent.Phenos)
# 
# #Plot each 'adaptive walk' on X and Y as a smooth regression, versus GEN.
# #Fig.3A Parallel
# Fig.3D <- ggplot(Fig.3D.Data, aes(x = 100*gen, y= TraitMean, colour = group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Phenotype") +
#   geom_smooth() + 
#   geom_hline(aes(yintercept = 0.25), colour = "black", linetype = "dashed") + 
#   geom_hline(aes(yintercept = -0.25), colour = "black", linetype = "dashed") + 
#   theme_ng1
# Fig.3D
# 
# #########
# ##FIG3E##
# #########
# 
# #Fig 3E genomic divergence in parallel with SGV+DNM or just DNM
# Fig3.Pos1.Muts.SGVDNM <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.251-0.251_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.Muts.SGVDNM <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.249-0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # Fig3.Neg1.Muts <- read.csv('data/parent_pop_K10000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.51--0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Pos1.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Pos1.Muts.SGVDNM)))
# MatrixList.Fig3.Pos2.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Pos2.Muts.SGVDNM)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3E.Pos1.Muts.SGVDNM <- lapply(MatrixList.Fig3.Pos1.Muts.SGVDNM, mutation.matrices) 
# CleanList.Fig3E.Pos2.Muts.SGVDNM <- lapply(MatrixList.Fig3.Pos2.Muts.SGVDNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3E.Parallel.Divergence.SGVDNM <- data.frame(mapply(between.pop.diversity.SGV, x = CleanList.Fig3.Pos1.Muts.SGVDNM, y = CleanList.Fig3.Pos1.Muts.SGVDNM, sgv.end = 145))
# 
# ## Fig 3E DNM 
# #Fig 3E genomic divergence in parallel with SGV+DNM or just DNM
# Fig3.Pos1.Muts.DNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.Muts.DNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.249-0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Pos1.Muts.DNM <- as.list(data.frame(t(Fig3.Pos1.Muts.DNM)))
# MatrixList.Fig3.Pos2.Muts.DNM <- as.list(data.frame(t(Fig3.Pos2.Muts.DNM)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3E.Pos1.Muts.DNM <- lapply(MatrixList.Fig3.Pos1.Muts.DNM, mutation.matrices) 
# CleanList.Fig3E.Pos2.Muts.DNM <- lapply(MatrixList.Fig3.Pos2.Muts.DNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3E.Parallel.Divergence.DNM <- data.frame(mapply(between.pop.diversity, x = CleanList.Fig3E.Pos1.Muts.DNM, y = CleanList.Fig3E.Pos2.Muts.DNM))
# 
# ## Make Fig 3E by plotting both manhattan distances...
# 
# Fig.3E.Data <- rbind(
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3E.Parallel.Divergence.SGVDNM$mapply.between.pop.diversity.SGV..x...CleanList.Fig3.Pos1.Muts.SGVDNM.., 
#              Group = rep("SGVDNM", 50)),
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3E.Parallel.Divergence.DNM$mapply.between.pop.diversity..x...CleanList.Fig3E.Pos1.Muts.DNM.., 
#              Group = rep("DNM", 50)))
# 
# Fig.3E <- ggplot(Fig.3E.Data, aes(x = 100*Gen, y= Dist, colour = Group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Genetic divergence") +
#   geom_smooth() + 
#   xlim(0,5000) +
#   ylim(0, 38) +
#   theme_ng1
# Fig.3E
# 
# #########
# ##FIG3F##
# #########
# 
# #Fig 3F SGVDNM
# Fig3.Neg1.Muts.SGVDNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt-0.249--0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Pos1.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Pos1.Muts.SGVDNM)))
# MatrixList.Fig3.Neg1.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Neg1.Muts.SGVDNM)))
# # matrix.list.Fig3.Neg1.Muts <- as.list(data.frame(t(Fig3.Neg1.Muts)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3.Pos1.Muts.SGVDNM <- lapply(MatrixList.Fig3.Pos1.Muts.SGVDNM, mutation.matrices) 
# CleanList.Fig3.Neg1.Muts.SGVDNM <- lapply(MatrixList.Fig3.Neg1.Muts.SGVDNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3F.Parallel.Divergence.SGVDNM <- data.frame(mapply(between.pop.diversity.SGV, x = CleanList.Fig3.Pos1.Muts.SGVDNM, y = CleanList.Fig3.Neg1.Muts.SGVDNM, sgv.end = 145))
# 
# ## Fig 3F DNM 
# #Fig 3F genomic divergence in divergent evolution with SGV+DNM or just DNM
# Fig3.Pos1.Muts.DNM <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Neg1.Muts.DNM <- fread('data/parent_pop_K2000_n2_B2_u0.001_alpha0.02_gens5000_founders2000_opt-0.249--0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Neg1.Muts.DNM <- as.list(data.frame(t(Fig3.Neg1.Muts.DNM)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3F.Neg1.Muts.DNM <- lapply(MatrixList.Fig3.Neg1.Muts.DNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3F.Parallel.Divergence.DNM <- data.frame(mapply(between.pop.diversity, x = CleanList.Fig3F.Pos1.Muts.DNM, y = CleanList.Fig3F.Neg1.Muts.DNM))
# 
# ## Make Fig 3F by plotting both manhattan distances...
# 
# Fig.3F.Data <- rbind(
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3F.Parallel.Divergence.SGVDNM$mapply.between.pop.diversity.SGV..x...CleanList.Fig3.Pos1.Muts.SGVDNM.., 
#              Group = rep("SGVDNM", 50)),
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3F.Parallel.Divergence.DNM$mapply.between.pop.diversity..x...CleanList.Fig3F.Pos1.Muts.DNM.., 
#              Group = rep("DNM", 50)))
# 
# # (force thru origin?)
# 
# Fig.3F <- ggplot(Fig.3F.Data, aes(x = 100*Gen, y= Dist, colour = Group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Genetic divergence") +
#   geom_smooth() + 
#   xlim(0,5000) +
#   ylim(0, 38) +
#   theme_ng1
# Fig.3F
# 
# ## Create multi-panel grid
# 
# Fig.3 <- plot_grid(Fig.3A, Fig.3B, Fig.3C, Fig.3D, Fig.3E, Fig.3F, ncol = 2, labels = c("a", "b", "c", "d", "e", "f"))
# 
# ############################
# ##Figure 4: Hybrid Fitness##
# ############################
# 
# Fig4.Parallel.SGVDNM.Hybrids <- read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_0.249-0.249_opt2_0.251-0.251_sgv.csv', header = F, col.names = c("X", "Y"))
# Fig4.Parallel.DNM.Hybrids <- read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_0.249-0.249_opt2_0.251-0.251_dnm.csv', header = F, col.names = c("X", "Y"))
# Fig4.Divergent.SGVDNM.Hybrids <- read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_-0.249--0.249_opt2_0.251-0.251_sgv.csv', header = F, col.names = c("X", "Y"))
# Fig4.Divergent.DNM.Hybrids <- read.csv('data/hybrid_phenos_K2000_n2_B2_u0.001_alpha0.02_gens5000_opt1_-0.249--0.249_opt2_0.251-0.251_dnm.csv', header = F, col.names = c("X", "Y"))
# 
# #Calculate hybrid fitness
# ## Parallel
# Fig4.Parallel.SGVDNM.Hybrid.Fitness <- data.frame(hybrid.fitness.n2.opt(Fig4.Parallel.SGVDNM.Hybrids$X, Fig4.Parallel.SGVDNM.Hybrids$Y, opt = 0.25, sigma.seln = 20))
# Fig4.Parallel.DNM.Hybrid.Fitness <- data.frame(hybrid.fitness.n2.opt(Fig4.Parallel.DNM.Hybrids$X, Fig4.Parallel.DNM.Hybrids$Y, opt = 0.25, sigma.seln = 20))
# 
# # Divergent
# 
# Fig4.Divergent.SGVDNM.Hybrid.Fitness <- data.frame(hybrid.fitness.n2.opt(Fig4.Divergent.SGVDNM.Hybrids$X, Fig4.Divergent.SGVDNM.Hybrids$Y, opt = 0.25, sigma.seln = 20))
# Fig4.Divergent.DNM.Hybrid.Fitness <- data.frame(hybrid.fitness.n2.opt(Fig4.Divergent.DNM.Hybrids$X, Fig4.Divergent.DNM.Hybrids$Y, opt = 0.25, sigma.seln = 20))
# 
# ## Make figure
# ## Fist thing first, make dataset
# 
# #Rename vars
# Fig4.P.SGV <- rename(Fig4.Parallel.SGVDNM.Hybrid.Fitness, Fitness = hybrid.fitness.n2.opt.Fig4.Parallel.SGVDNM.Hybrids.X..Fig4.Parallel.SGVDNM.Hybrids.Y..)
# Fig4.P.DNM <- rename(Fig4.Parallel.DNM.Hybrid.Fitness, Fitness = hybrid.fitness.n2.opt.Fig4.Parallel.DNM.Hybrids.X..Fig4.Parallel.DNM.Hybrids.Y..)
# Fig4.D.SGV <- rename(Fig4.Divergent.SGVDNM.Hybrid.Fitness, Fitness = hybrid.fitness.n2.opt.Fig4.Divergent.SGVDNM.Hybrids.X..Fig4.Divergent.SGVDNM.Hybrids.Y..)
# Fig4.D.DNM <- rename(Fig4.Divergent.DNM.Hybrid.Fitness, Fitness = hybrid.fitness.n2.opt.Fig4.Divergent.DNM.Hybrids.X..Fig4.Divergent.DNM.Hybrids.Y..)
# 
# Fig4.P.SGV$Class <- rep("SGV + DNM", nrow(Fig4.P.SGV))
# Fig4.P.DNM$Class <- rep("DNM", nrow(Fig4.P.DNM))
# Fig4.D.SGV$Class <- rep("SGV + DNM", nrow(Fig4.D.SGV))
# Fig4.D.DNM$Class <- rep("DNM", nrow(Fig4.D.DNM))
# 
# Fig4.P.SGV$ParDiv <- rep("Par", nrow(Fig4.P.SGV))
# Fig4.P.DNM$ParDiv <- rep("Par", nrow(Fig4.P.DNM))
# Fig4.D.SGV$ParDiv <- rep("Div", nrow(Fig4.D.SGV))
# Fig4.D.DNM$ParDiv <- rep("Div", nrow(Fig4.D.DNM))
# 
# 
# Fig4.Data <- rbind(Fig4.P.SGV, Fig4.P.DNM, Fig4.D.SGV, Fig4.D.DNM)
# Fig4.Data$ParDiv <- as.factor(Fig4.Data$ParDiv)
# Fig4.Data$Class <- as.factor(Fig4.Data$Class)
# 
# 
# Fig4.Data.Plot <- Rmisc::summarySE(Fig4.Data, measurevar = "Fitness", groupvars = c("Class", "ParDiv"))
# Fig4.Data.Plot$ParDiv <- relevel(Fig4.Data.Plot$ParDiv, ref = "Par")
# Fig4.Data.Plot$Class <- relevel(Fig4.Data.Plot$Class, ref = "DNM")
# 
# 
# Fig4 <- ggplot(Fig4.Data.Plot, aes(group = ParDiv, x= factor(Class, labels=c("DNM", "SGV")), y=Fitness)) +
#   labs(x = NULL, y = ("Fitness")) +
#   geom_errorbar(aes(ymin = Fitness - sd, ymax = Fitness + sd), width = 0.3, position = position_dodge(0.3)) +
#   geom_point(position = position_dodge(0.3), colour = "black", stat = "identity", aes(shape = ParDiv, size = 4)) +
#   scale_fill_manual(values=c("white", "black")) +
#   scale_shape_manual(values=c(24, 19)) +
#   geom_line(position = position_dodge(0.3), aes(group = ParDiv)) + 
#   theme_ng1
# Fig4
# 
# ggsave(filename = '/Users/Ken/Dropbox/!Ph.D./!SVS_Standing-Variation-Speciation/SVS_Poster/SVS_poster_figs/Fig2.pdf' , plot = Fig4)
#   theme_ng1
# Fig.3D
# 
# #########
# ##FIG3E##
# #########
# 
# #Fig 3E genomic divergence in parallel with SGV+DNM or just DNM
# Fig3.Pos1.Muts.SGVDNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.251-0.251_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.Muts.SGVDNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.249-0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# # Fig3.Neg1.Muts <- read.csv('data/parent_pop_K10000_n2_B2_u0.001_alpha0.02_gens5000_opt-0.51--0.51.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Pos1.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Pos1.Muts.SGVDNM)))
# MatrixList.Fig3.Pos2.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Pos2.Muts.SGVDNM)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3E.Pos1.Muts.SGVDNM <- lapply(MatrixList.Fig3.Pos1.Muts.SGVDNM, mutation.matrices) 
# CleanList.Fig3E.Pos2.Muts.SGVDNM <- lapply(MatrixList.Fig3.Pos2.Muts.SGVDNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3E.Parallel.Divergence.SGVDNM <- data.frame(mapply(between.pop.diversity.SGV, x = CleanList.Fig3.Pos1.Muts.SGVDNM, y = CleanList.Fig3.Pos1.Muts.SGVDNM, sgv.end = 145))
# 
# ## Fig 3E DNM 
# #Fig 3E genomic divergence in parallel with SGV+DNM or just DNM
# Fig3.Pos1.Muts.DNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Pos2.Muts.DNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.249-0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Pos1.Muts.DNM <- as.list(data.frame(t(Fig3.Pos1.Muts.DNM)))
# MatrixList.Fig3.Pos2.Muts.DNM <- as.list(data.frame(t(Fig3.Pos2.Muts.DNM)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3E.Pos1.Muts.DNM <- lapply(MatrixList.Fig3.Pos1.Muts.DNM, mutation.matrices) 
# CleanList.Fig3E.Pos2.Muts.DNM <- lapply(MatrixList.Fig3.Pos2.Muts.DNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3E.Parallel.Divergence.DNM <- data.frame(mapply(between.pop.diversity, x = CleanList.Fig3E.Pos1.Muts.DNM, y = CleanList.Fig3E.Pos2.Muts.DNM))
# 
# ## Make Fig 3E by plotting both manhattan distances...
# 
# Fig.3E.Data <- rbind(
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3E.Parallel.Divergence.SGVDNM$mapply.between.pop.diversity.SGV..x...CleanList.Fig3.Pos1.Muts.SGVDNM.., 
#              Group = rep("SGVDNM", 50)),
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3E.Parallel.Divergence.DNM$mapply.between.pop.diversity..x...CleanList.Fig3E.Pos1.Muts.DNM.., 
#              Group = rep("DNM", 50)))
# 
# Fig.3E <- ggplot(Fig.3E.Data, aes(x = 100*Gen, y= Dist, colour = Group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Genetic divergence") +
#   geom_smooth() + 
#   xlim(0,5000) +
#   ylim(0, 38) +
#   theme_ng1
# Fig.3E
# 
# #########
# ##FIG3F##
# #########
# 
# #Fig 3F SGVDNM
# Fig3.Neg1.Muts.SGVDNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt-0.249--0.249_sgv.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Pos1.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Pos1.Muts.SGVDNM)))
# MatrixList.Fig3.Neg1.Muts.SGVDNM <- as.list(data.frame(t(Fig3.Neg1.Muts.SGVDNM)))
# # matrix.list.Fig3.Neg1.Muts <- as.list(data.frame(t(Fig3.Neg1.Muts)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3.Pos1.Muts.SGVDNM <- lapply(MatrixList.Fig3.Pos1.Muts.SGVDNM, mutation.matrices) 
# CleanList.Fig3.Neg1.Muts.SGVDNM <- lapply(MatrixList.Fig3.Neg1.Muts.SGVDNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3F.Parallel.Divergence.SGVDNM <- data.frame(mapply(between.pop.diversity.SGV, x = CleanList.Fig3.Pos1.Muts.SGVDNM, y = CleanList.Fig3.Neg1.Muts.SGVDNM, sgv.end = 145))
# 
# ## Fig 3F DNM 
# #Fig 3F genomic divergence in divergent evolution with SGV+DNM or just DNM
# Fig3.Pos1.Muts.DNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt0.251-0.251_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# Fig3.Neg1.Muts.DNM <- fread('data/parent_pop_K1000_n2_B2_u0.001_alpha0.02_gens5000_founders1000_opt-0.249--0.249_dnm.csv', header = F, check.names = F, na.strings = c("[", "]"))
# 
# # Make each generation an item in a list
# MatrixList.Fig3.Neg1.Muts.DNM <- as.list(data.frame(t(Fig3.Neg1.Muts.DNM)))
# 
# #Use 'mutation.matrices' function to make each item in a list a matrix that can be used to calculate Euclidean distance
# CleanList.Fig3F.Neg1.Muts.DNM <- lapply(MatrixList.Fig3.Neg1.Muts.DNM, mutation.matrices) 
# 
# # Determine between population diversity (takes a few minutes)
# Fig.3F.Parallel.Divergence.DNM <- data.frame(mapply(between.pop.diversity, x = CleanList.Fig3F.Pos1.Muts.DNM, y = CleanList.Fig3F.Neg1.Muts.DNM))
# 
# ## Make Fig 3F by plotting both manhattan distances...
# 
# Fig.3F.Data <- rbind(
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3F.Parallel.Divergence.SGVDNM$mapply.between.pop.diversity.SGV..x...CleanList.Fig3.Pos1.Muts.SGVDNM.., 
#              Group = rep("SGVDNM", 50)),
#   data.frame(Gen = c(1:50), 
#              Dist = Fig.3F.Parallel.Divergence.DNM$mapply.between.pop.diversity..x...CleanList.Fig3F.Pos1.Muts.DNM.., 
#              Group = rep("DNM", 50)))
# 
# # (force thru origin?)
# 
# Fig.3F <- ggplot(Fig.3F.Data, aes(x = 100*Gen, y= Dist, colour = Group)) +
#   geom_point() + 
#   labs(x = "Generation #",
#        y = "Genetic divergence") +
#   geom_smooth() + 
#   xlim(0,5000) +
#   ylim(0, 38) +
#   theme_ng1
# Fig.3F
# 
# ## Create multi-panel grid
# 
# Fig.3 <- plot_grid(Fig.3A, Fig.3B, Fig.3C, Fig.3D, Fig.3E, Fig.3F, ncol = 2, labels = c("a", "b", "c", "d", "e", "f"))
# 
# ############################
# ##Figure 4: Hybrid Fitness##
# ############################
# 
# ##Currently not as big of a difference as I'd like between SGV / DNM sims in parallel... it's probably because there is lower parallelism than expected...
# 
# Fig4.Parallel.SGVDNM.Hybrids <- read.csv('data/hybrid_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt1_0.249-0.249_opt2_0.251-0.251_sgv.csv', header = F, col.names = c("X", "Y"))
# Fig4.Parallel.DNM.Hybrids <- read.csv('data/hybrid_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt1_0.249-0.249_opt2_0.251-0.251_dnm.csv', header = F, col.names = c("X", "Y"))
# Fig4.Divergent.SGVDNM.Hybrids <- read.csv('data/hybrid_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt1_-0.249--0.249_opt2_0.251-0.251_sgv.csv', header = F, col.names = c("X", "Y"))
# Fig4.Divergent.DNM.Hybrids <- read.csv('data/hybrid_phenos_K1000_n2_B2_u0.001_alpha0.02_gens5000_opt1_-0.249--0.249_opt2_0.251-0.251_dnm.csv', header = F, col.names = c("X", "Y"))
# 
# #Calculate hybrid fitness
# Fig4.Parallel.SGVDNM.Hybrid.Fitness <- hybrid.fitness.n2.opt(Fig4.Parallel.SGVDNM.Hybrids$X, Fig4.Parallel.SGVDNM.Hybrids$Y, opt = 0.25)
# Fig4.Parallel.DNM.Hybrid.Fitness <- hybrid.fitness.n2.opt(Fig4.Parallel.DNM.Hybrids$X, Fig4.Parallel.DNM.Hybrids$Y, opt = 0.25)
# 
# #Parallel hybrids
# 
# plot(Fig4.Parallel.SGVDNM.Hybrids$X, Fig4.Parallel.SGVDNM.Hybrids$Y, col = "red")
# points(Fig4.Parallel.DNM.Hybrids$X, Fig4.Parallel.DNM.Hybrids$Y, col = "green")
# 
# sd(Fig4.Parallel.DNM.Hybrids$Y)
# 
# ######################################################################
# ##SUPPLEMENTARY FIGURES##
# ######################################################################
# 
# #FigSX -- Adaptation to parallel and divergent optima with ONLY SGV, and 'genomic' divergence.