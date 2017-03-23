## Standing Genetic Variation & Speciation R Analysis script
## Author Ken A. Thompson

#load packages
library(tidyverse)


## Testing script for plotting adaptive walk
#Generate vectors

mut.d1 <- rnorm(10) %>% 
  as.vector() %>% 
  abs() %>% 
  sort()
  
mut.d2 <- rnorm(10) %>% 
  as.vector() %>% 
  abs() %>% 
  sort()


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
  geom_segment(aes(xend=c(tail(mut.d1, n=-1), NA), yend=c(tail(mut.d2, n=-1), NA)),
               arrow=arrow(length=unit(0.3,"cm"), type = "open")) + 
  theme_ng1
adapt.walk


