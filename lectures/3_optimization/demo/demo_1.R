library(MASS)
library(rgl)
library(R.matlab)

source('kernels.R')
source('models.R')
source('likelihood.R')


# get data
data <- readMat('data_nonoise.mat')[[1]]
X <- data[,1:2]
F <- data[,4]

xmin <- min(X[,1])
xmax <- max(X[,1])

ymin <- min(X[,2])
ymax <- max(X[,2])

xgrid <- seq(xmin, xmax, length=30)
ygrid <- seq(ymin, ymax, length=30)

Xgrid <- expand.grid(xgrid,ygrid)



##################################
## Multivariate GPR

## look at the data
plot3d(X[,1], X[,2], F)

## what prior seem appropriate ?
open3d()
K <- kExp(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

open3d()
K <- kMat32(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

open3d()
K <- kMat52(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')

open3d()
K <- kGauss(Xgrid,Xgrid,param=c(1,3000,3000))
Z <- mvrnorm(1,rep(0,nrow(K)),K)
persp3d(xgrid,ygrid,Z,col='wheat')


####################################
## Parameter estimation

## plot log likelihood
sig2_grid <- seq(1e-4, 1e-3, length=30)
theta_grid <- seq(500, 3000, length=30)

PARAM <- as.matrix(expand.grid(sig2_grid, theta_grid))
LL <- rep(0, nrow(PARAM))

for(i in 1:nrow(PARAM)){
    LL[i] <- logLikelihood(c(PARAM[i,],PARAM[i,2]),kGauss,X,F)
}

persp3d(sig2_grid, theta_grid, LL, col='red')

## optimisation
indmax <- which.max(LL)
paropt <- c(PARAM[indmax,], PARAM[indmax,2])
spheres3d(paropt[1], paropt[2], LL[indmax], radius=20)

## model associated with optimal parameters
pred <- predGPR(Xgrid, X, F,
                kern=kGauss, param=paropt)

pred_sd <- pmax(rep(0,ncol(pred$cov)), diag(pred$cov))
pred_upper95 <- pred$mean + 2*sqrt(pred_sd)
pred_lower95 <- pred$mean - 2*sqrt(pred_sd)

persp3d(xgrid, ygrid, pred$mean, col='red')
surface3d(xgrid, ygrid, pred_lower95, col='wheat', alpha=.8)
surface3d(xgrid, ygrid, pred_upper95, col='wheat', alpha=.8)
spheres3d(X[,1], X[,2], F, radius=200)

##################################
## Multivariate GPR with noise

Kn <- kExp(X, X, param=c(1e-4, 2000, 2000)) + kWhite(X, X, param=1e-5)
Fn <- mvrnorm(1,F,Kn)

## look at the data
plot3d(X[,1], X[,2], Fn)

## construct appropriate kernel for the signal
kSignal <- function(x, y, param){
    sig2 <- abs(param[1])
    theta <- rep(param[2], 2)
    decay <- rep(param[3], 2)

    loc <- matrix(c(366000, 7650000), ncol=2)
    locx <-exp(-.5*mydist(x, loc, decay))
    locy <-exp(-.5*mydist(y, loc, decay))
    Locx_mat <- matrix(locx, nrow=nrow(x), ncol=nrow(y))
    Locy_mat <- t(matrix(locy, nrow=nrow(y), ncol=nrow(x)))

    kgauss <- exp(-.5*mydist(x, y, theta)^2)
 
    kern <- sig2 * Locx_mat * Locy_mat * kgauss

    return(kern)
}

KS <- kSignal(Xgrid, Xgrid, c(1e-3, 4000, 3000))
Z <- mvrnorm(1, rep(0,nrow(KS)), KS)
persp3d(xgrid, ygrid, Z, col='wheat')

## construct appropriate kernel for the noise
kNoise <- function(x, y, param){
    kW <- kWhite(x, y, param[1])
    kE <- kExp(x, y, param[c(2,3,3)])
    return(kW + kE)
}

KN <- kNoise(Xgrid, Xgrid, c(1e-5, 1e-4, 3000))
ZN <- mvrnorm(1, rep(0,nrow(KN)), KN)
persp3d(xgrid, ygrid, ZN, col='wheat')

## Optimize the log-likelihood
params <- c(1e-3, 4000, 6000, 1e-5, 1e-4, 2000)
logLikelihood(params, kSignal, X, Fn, kNoise, num_param_noise=3)

opt <- optim(params, logLikelihood, control=list(fnscale=-1), 
             kern=kSignal, Xd=X, F=Fn, kernNoise=kNoise, num_param_noise=3)

## build a GPR model
pred <- predGPR(Xgrid, X, Fn,
                kern=kSignal, param=opt$par[1:3],
                kernNoise=kNoise, paramNoise=opt$par[4:6])

pred_sd <- pmax(rep(0,ncol(pred$cov)), diag(pred$cov))
pred_upper95 <- pred$mean + 2*sqrt(pred_sd)
pred_lower95 <- pred$mean - 2*sqrt(pred_sd)

persp3d(xgrid, ygrid, pred$mean, col='red')
surface3d(xgrid, ygrid, pred_lower95, col='wheat', alpha=.5)
surface3d(xgrid, ygrid, pred_upper95, col='wheat', alpha=.5)
spheres3d(X[,1], X[,2], Fn, radius=200)
