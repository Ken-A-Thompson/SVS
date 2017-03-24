## Standing Genetic Variation & Speciation R Analysis script
## Author Ken A. Thompson

#load packages
#install.packages('tidyverse')
library(tidyverse)


## Testing script for plotting adaptive walk
#Generate vectors

mut.d1 <- c(0, 
            rnorm(19) %>% 
  as.vector() %>% 
  abs() %>% 
  sort())
  
mut.d2 <- c(0, 
            rnorm(19) %>% 
  as.vector() %>% 
  abs() %>% 
  sort())

#Generate a 'hybrid'

## Create 'parents'
anc.d1 <- rep(0, length(mut.d1))
anc.d2 <- rep(0, length(mut.d2))

#Generate Hybrids
muts.d1 <- diff(mut.d1, lag = 1)
muts.d2 <- diff(mut.d2, lag = 1)

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
adapt.walk <- ggplot(df, aes(x = mut.d1, y= mut.d2)) +
  geom_point() + 
  labs(x = "Trait 1", y = "Trait 2") +
  geom_segment(aes(xend=c(tail(mut.d1, n=-1), NA), yend=c(tail(mut.d2, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open")) + 
  theme_ng1
adapt.walk

#Plot hybrids
## Create hybrid dataframe

hyb.df <- data.frame(hyb.d1, hyb.d2)

adapt.walk + geom_point(data = hyb.df, aes(x = hyb.d1, y = hyb.d2), alpha = 0.25)

