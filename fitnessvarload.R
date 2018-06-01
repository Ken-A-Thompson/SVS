install.packages("R2Cuba")
library("R2Cuba")

install.packages("mvtnorm")
library("mvtnorm")

integrand <- function(arg, theta, lambda) {
  x <- arg[1] #trait 1
  y <- arg[2] #trait 2
  x1 <- 1 #optimum x in enviro 1
  y1 <- 0 #optimum y in enviro 1
  x2 <- cos(theta*pi/180) #optimum x in enviro 2
  y2 <- sin(theta*pi/180) #optimum x in enviro 2
  pxy <- dmvnorm( arg, mean = c((x1+x2)/2, (y1+y2)/2), sigma = lambda*diag(2), log=FALSE ) #distribution of traits (normal with mean between optima and variance lambda)
  wxy1 <- 2*pi*dmvnorm( arg, mean = c(x1, y2), sigma = diag(2), log=FALSE ) #fitness in enviro 1
  wxy2 <- 2*pi*dmvnorm( arg, mean = c(x1, y2), sigma = diag(2), log=FALSE ) #fitness in enviro 2
  novar <- 2*pi*dmvnorm( c((x1+x2)/2, (y1+y2)/2), mean = c(x1, y2), sigma = diag(2), log=FALSE ) #mean fitness if no variation
  ff <- max(wxy1,wxy2) * pxy / novar; #max fitness times density divided by no variance fitness (effect of variance on mean fitness)
  return(ff)
} # end integrand

points <- function(lambda, theta){
  varload <- -log(1/(1+lambda)) #variance load
  deltaw <- cuhre(ndim=2, ncomp=1, integrand, theta=theta, lambda=lambda, lower=rep(-100,2), upper=rep(100,2), rel.tol=1e-3, abs.tol=1e-12, flags=list(verbose=2)) #effect on mean fitness
  return(c(theta,varload,deltaw$value)) #return angle, variance load, and effect of variance on mean fitness
}

points(0.001,60)
