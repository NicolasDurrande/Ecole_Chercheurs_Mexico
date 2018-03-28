library(DiceKriging)
library(DiceView)
library(DiceOptim)
source("utilities_volcan.R")

#### Question 1 ####


#### Question 2 ####
m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", nugget=1e-8)

#### Question 3 ####
print(m)

#### Question 4 ####


#### Question 5 ####
err_std <- ... # calcul des résidus standardisés

# test des intervalles de confiance
x <- seq(-3, 3, 0.03)
hist(err_std, freq=FALSE, xlim=c(min(err_std, -3), max(err_std, 3)))
lines(x, dnorm(x))

#####################
## optimisation globale
#####################

#### Question 6 ####


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

#### Question 10 ####


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
