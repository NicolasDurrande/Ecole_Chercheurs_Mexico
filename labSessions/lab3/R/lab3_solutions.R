library(DiceKriging)
library(DiceView)
library(DiceOptim)

#### Question 1 ####
load('XY_volcano.Rdata')

# test des données chargées
par(mfrow=c(1,5))
for(i in 1:5){
  plot(X[,i], Y)
}
par(mfrow=c(1,1))

#### Question 2 ####
m1 <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", nugget=1e-8)

m2 <- km(y ~ .^2 + x.1^2 + x.2^2 + x.3^2 + x.4^2 + x.5^2, design=data.frame(x=X), response=data.frame(y=Y), 
         covtype="gauss", nugget=1e-8, lower=rep(0.05, 5))

#m2 <- km(y ~ 1, design=data.frame(x=X), response=data.frame(y=Y), 
#         covtype="matern5_2", nugget=1e-8, lower=rep(0.05, 5))

m3 <- km(y~., design=data.frame(x=X), response=data.frame(y=Y), 
         covtype="matern5_2", nugget=1e-8)

# calculer la (log) vraisemblance
c(logLik(m1), logLik(m2), logLik(m3))

#### Question 3 ####
print(m2)

#### Question 4 ####
Q2 <- function(Y, Ypred){
  1 - sum((Y - Ypred)^2) / sum((Y - mean(Y))^2)
}

loo1 <-leaveOneOut.km(m1, type='UK')
loo2 <-leaveOneOut.km(m2, type='UK')
loo3 <-leaveOneOut.km(m3, type='UK')

c(Q2(Y, loo1$mean), Q2(Y, loo2$mean), Q2(Y, loo3$mean))

#### Question 5 ####
res1_std <- (Y-loo1$mean) / loo1$sd
res2_std <- (Y-loo2$mean) / loo2$sd
res3_std <- (Y-loo3$mean) / loo3$sd

# test des intervalles de confiance
x <- seq(-3, 3, 0.03)

par(mfrow=c(1,3))
for(res in list(res1_std, res2_std, res3_std)){
  hist(res, freq=FALSE, xlim=c(min(res, -3), max(res, 3)))
  lines(x, dnorm(x))
}
par(mfrow=c(1,1))

m <- m2

#####################
## optimisation globale
source("utilities_volcan.R")

#### Question 6 ####
res <- EGO.nsteps(model=m, fun=compute_wls, nsteps=50, lower=rep(0, 5),
                  upper=rep(1,5), control=list(print.level=0))

print(res$lastmodel)
ybest <- min(res$lastmodel@y)

#### Question 7 ####
n <- length(Y)
p <- length(res$lastmodel@y)
best <- cummin(res$lastmodel@y)

plot(1:p, res$lastmodel@y)
lines(1:p, best)
abline(v=n+.5, lty=2)

#### Question 8 ####

cols=c(rep("black", n), rep("red", p))
cols[which(res$lastmodel@y < -1.6)] <- "blue"
X2 <- res$lastmodel@X
pairs(X2, col=cols)

#### Question 9 ####
xbest <- res$lastmodel@X[which.min(res$lastmodel@y),]
sectionview(res$lastmodel, center = xbest, ylim=c(-2.7, 2))

#### Question 10 ####
mc <- res$lastmodel
resc <- EGO.nsteps(model=mc, fun=compute_wls, nsteps=50, lower=rep(0, 5),
                  upper=rep(1,5), control=list(print.level=0))

print(resc$lastmodel)
ybest <- min(resc$lastmodel@y)

#### Question 11 ####
library(lhs)
getmean <- function(newdata, m) {
  pred <- predict(object=m, newdata=newdata, type="UK")
  return(pred$mean)
}
X1 <- data.frame(randomLHS(10000,5))
X2 <- data.frame(randomLHS(10000,5))
colnames(X1) <- colnames(X2) <- colnames(m@X)
res2 <- soboljansen(model = getmean, X1=X1, X2=X2, nboot = 50, conf = 0.95, m=m)
plot(res2)

X1 <- data.frame(randomLHS(1000,5))
X2 <- data.frame(randomLHS(1000,5))
candidate <- data.frame(randomLHS(100,5))
colnames(X1) <- colnames(X2) <- colnames(candidate) <- colnames(m@X)
res <- sobolGP(model = m, type="UK", MCmethod="soboljansen",
               X1=X1, X2=X2, nsim = 20, nboot=50, sequential = TRUE, candidate=candidate)

plot(res)
