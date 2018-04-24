#setwd("/Users/mmosmond/Documents/PHD/SVS")
data <- read.csv("data/burnins/PlotBurn_n2_K10000_alpha0.1_B2_u0.0010_sigma0.1.csv", header=FALSE) #read csv

library(ggplot2)
ggplot(data=data, aes(x=V2, y=V3)) + 
  geom_line(aes(group=V1)) + geom_point()
