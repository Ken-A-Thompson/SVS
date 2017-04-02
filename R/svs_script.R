## Standing Genetic Variation & Speciation R Analysis script
## Author Ken A. Thompson
### 02 Apr 2017

#load packages
#install.packages('tidyverse')
#install.packages('splitstackshape')
library(splitstackshape)
library(tidyverse)

## Testing script for plotting adaptive walk
#Generate vectors

mut.d1 <- c(0, 
            rnorm(20)) %>%
  abs() %>% 
  sort() %>% 
  as.vector()
  
mut.d2 <- c(0, 
            rnorm(20)) %>%
  abs() %>% 
  sort() %>% 
  as.vector()

# Generate a 'hybrid'

## Create 'parents'
anc.d1 <- rep(0, length(mut.d1))
anc.d2 <- rep(0, length(mut.d2))

#Generate Hybrids
muts.d1 <- rev(diff(mut.d1, lag = 1))
muts.d2 <- rev(diff(mut.d2, lag = 1))

muts.d1.vec <- c(0, muts.d1)
muts.d2.vec <- c(0, muts.d2)

hyb.d1 <- replicate(1000, sum(sample(c(muts.d1.vec, anc.d1), length(muts.d1.vec),replace=F))) #Sum of 'trait value' for D1
hyb.d2 <- replicate(1000, sum(sample(c(muts.d2.vec, anc.d2), length(muts.d2.vec),replace=F))) #Sum of 'trait value' for D2


## Generate figure for 'adaptive' walk
### Run theme

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

### End theme 

df <- data.frame(mut.d1, mut.d2) # Generate data frame for ggplot with both mutation vectors

### plot it 
adapt.walk <- ggplot(df, aes(x = rev(mut.d1), y= rev(mut.d2))) +
  geom_point() + 
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(rev(mut.d1), n=-1), NA), yend=c(tail(rev(mut.d2), n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open")) + 
  theme_ng1
adapt.walk

#Plot hybrids
## Create hybrid dataframe

hyb.df <- data.frame(hyb.d1, hyb.d2)

adapt.walk + geom_point(data = hyb.df, aes(x = hyb.d1, y = hyb.d2), alpha = 0.25)


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








