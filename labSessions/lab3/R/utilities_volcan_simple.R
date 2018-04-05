###########################
# Utility functions for
# Inversion of a punctual displacements source from 3D data
# The data used are under-sampled with the quadtree method (irregular grid)
#
# Rodolphe Le Riche, Nicolas Durrande, Valerie Cayol, Victor Picheny
#
#
# Notes: 
# * use of global variables, starting with Glb_...
###########################

####### load utilities ##########################

# rm(list=ls()) #  cleaning up

library(R.matlab)
source("./mogi_3D.R")
source("./wls_ulos.R")
source('kernels.R')

###### input for variables identification ###########
n2i <- list(xs=1, ys=2, zs=3, a=4, p=5) # name to index for variables
varnames <- c("xs","ys","zs","a","p")
nbvar <- 5

# optimum
xstar <- NA
xstar[n2i$xs] <- 367000 # X location of source in m UTM coordinates
xstar[n2i$ys] <- 7650300 # Y location of source in m UTM
xstar[n2i$zs] <- 0 # Elevation of source with respect to sea level in m
xstar[n2i$a] <- 500 # source radius in m
xstar[n2i$p] <- 20 # Source overpressure in MPa
xstar_norm <- c(0.75, 0.65, 0.75, 0.4736842, 0.52)

# order of magnitude of the variables 
# (different from 0, to be used when scaling, setting bounds...)
xmag <- NA
xmag[n2i$xs] <- 367000 # X location of source in UTM coordinates
xmag[n2i$ys] <- 7650300 # Y location of source in UTM
xmag[n2i$zs] <- 1000 # Elevation of source with respect to sea level
xmag[n2i$a] <- 500 # source radius
xmag[n2i$p] <- 100 # Source overpressure in MPa

# bounds on variables
xmax <- NA
xmin <- NA
xmin[n2i$xs]<-366000
xmax[n2i$xs]<-368000
xmin[n2i$ys]<-7650000
xmax[n2i$ys]<-7651000
xmin[n2i$zs]<- -3000
xmax[n2i$zs]<-1000
xmin[n2i$a]<- 100
xmax[n2i$a]<- 1000
xmin[n2i$p]<- -500
xmax[n2i$p]<-500

Glb_var <<- list(n2i=n2i,nbvar=nbvar,xmag=xmag,xmax=xmax,xmin=xmin) # always useful stuff

####### load data ###########

data <- readMat('data_nonoise.mat') # TODO: do a read from cvs version (rodo)
Glb_xi <<- as.matrix(data$locdata[,1])
Glb_yi <<- as.matrix(data$locdata[,2])
Glb_zi <<- as.matrix(data$locdata[,3])
Glb_ulos <<- as.matrix(data$locdata[,4])

# calculate data Covariance matrix, store it in a Global variable
# covariance from exponential kernel, var = 5e-4m2, cor_length = 850 m
# and invert it
Xdata <- data$locdata[,1:2] # z's are not accounted for in Xdata

Glb_CXinv <<- solve(kExp(Xdata,Xdata,c(5e-4,850,850))) # calculated once for all, used in wls_ulos
rm(data)
rm(Xdata)

#######  useful functions #############################
# Scale from [0 1] to [xmin xmax]
unnorm_var <- function(Xnorm){
  if (is.null(dim(Xnorm))) Xnorm <- matrix(data = Xnorm, nrow=1) # numeric vector
  nbrep <- nrow(Xnorm)
  Xu <- matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar) + 
    Xnorm * matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)
  colnames(Xu) <- varnames
  return(Xu)
}

# Scale from [xmin xmax] to [0 1] 
norm_var <- function(X){
  if (is.null(dim(X))) X <- matrix(data = X, nrow=1)
  nbrep <- nrow(X)
  Xn <- (X - matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar)) / 
    matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)  
  colnames(Xn) <- varnames
  return(Xn)
}

# normalize the output so that it is centered with a unit std dev
# because wls ranges from 0 to 10^9, do a scaling in log(1+wls)
# wls normalization
normalizeWLS <- function(awls){
  lawls <- log(1+awls)
  return((lawls - 8.54)/3.2)
}

# objective function
compute_wls <- function(x) {
  if (is.null(dim(x))) {
    x <- matrix(data = x, nrow=1)
  }  else if (ncol(x)!=nbvar) {
    x <- t(x)
  }
  xx <- unnorm_var(Xnorm=x)
  y <- apply(xx, 1, wls_ulos)
  return(normalizeWLS(y))
}
