#fraction of beneficial mutations not beneficial in the other environment, at the beginning of the walk, for a number of dimensions (m), plotted over angles between optima  
m <- 2
curve(1 - pbeta(cos(x*(pi/180)/2)^2, (1+m)/2, 1/2), 0, 180, xlab="angle (degrees)", ylab="fraction non-overlap")
m <- 4
curve(1 - pbeta(cos(x*(pi/180)/2)^2, (1+m)/2, 1/2), 0, 180, add=T)
m <- 8
curve(1 - pbeta(cos(x*(pi/180)/2)^2, (1+m)/2, 1/2), 0, 180, add=T)

#fraction of beneficial mutations also beneficial in the other environment, plotted across the walk, for a given dimension (m) and a few angles (theta)
m <-6
theta<-0
curve(pbeta(1 + (cos(theta*(pi/180)) - 1)/(2*(1-x)^2), (1+m)/2, 1/2), 0, 1, xlab="fraction traversed", ylab="fraction overlap", ylim=c(0,1))
theta<-30
curve(pbeta(1 + (cos(theta*(pi/180)) - 1)/(2*(1-x)^2), (1+m)/2, 1/2), 0, 1, add=T)
theta<-60
curve(pbeta(1 + (cos(theta*(pi/180)) - 1)/(2*(1-x)^2), (1+m)/2, 1/2), 0, 1, add=T)
