# set the working directory such that it contains the volcano file
# setwd('/home/nicolas/Documents/courses/VolcanTestCase') 
setwd('/home/vpicheny/Documents/Ecole chercheur/VolcanTestCase') 
# load the libraries we will be using
library(DiceDesign)
source('utilities_volcan.R')

# set some variables
n <- 100           # nb of points in the DoE
d <- 5             # dimension of the input space

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

###### Q1 ######

unif_doe <- function(n, d){
    # inputs:
    #    n (int): number of points in the DoE
    #    d (int): dimension of the input space
    # outputs:
    #    a n*d matrix of points uniformly sampled in [0,1]^d
    X <- matrix(runif(n*d), ncol=d)
    return(X)
}

X_rand <- unif_doe(n, d)
pairs(X_rand, diag.panel=panel.hist)


###### Q2 ######

lhs_doe <- function(n, d){
    # inputs:
    #    n (int): number of points in the DoE
    #    d (int): dimension of the input space
    # outputs:
    #    a n*d matrix of a Latin Hypercube desing in [0,1]^d
    X <- matrix(0, nrow=n, ncol=d)
    x <- (1:n)/n-1/(2*n)
    for(i in 1:d){
        X[,i] <- sample(x)
    }
    return(X)
}

X_lhs <- lhs_doe(n, d)
pairs(X_lhs, diag.panel=panel.hist)
summary(X_lhs)


###### Q3 ######
library(DiceDesign)
X <- lhsDesign(n, d)$design
colnames(X) <- 1:5
Xopt <- maximinESE_LHS(design=X, inner_it=10, it=1)
plot(Xopt$critValues,type="l")
X_lhsopt <- Xopt$design
pairs(X_lhsopt,  diag.panel=panel.hist)


###### Q4 ######

X_faure <- runif.faure(n, d)$design
pairs(X_faure, diag.panel=panel.hist)


###### Q5 ######

X <- X_rand
hist(X)
plot(X[,1], rep(1, n), type="h")
pairs(X)

Xlist <- list(X_rand, X_lhs, X_lhsopt, X_faure)

unlist(lapply(Xlist, mindist))
unlist(lapply(Xlist, meshRatio))

mindist(X)
meshRatio(X)

###### Q6 ######

X <- X_lhsopt

Y <- compute_wls(X)


###### Q6 ######

par(mfrow=c(1,5))
for(i in 1:5){
    plot(X[,i], Y, xlab=paste0("X",i), ylab="y")
}
par(mfrow=c(1,1))
