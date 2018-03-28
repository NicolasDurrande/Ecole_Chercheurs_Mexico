library(DiceKriging)
library(DiceView)

#### Question 1 ####
load('XY_volcano.Rdata')

# test des données chargées
par(mfrow=c(1,5))
for(i in 1:5){
  plot(X[,i], Y)
}
par(mfrow=c(1,1))

#### Question 2 ####
# faire du multistart !
m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", nugget=1e-8)

print(m)

# calculer la (log) vraisemblance
logLik(m)

#### Question 2 ####


#### Question 9 ####
res <-leaveOneOut.km(m, type='SK')
res_std <- (Y-res$mean) / res$sd

# test de la moyenne
Q2 <- 1 - sum((Y - res$mean)^2) / sum((Y - mean(Y))^2)

# test des intervalles de confiance
hist(res_std, freq=FALSE, xlim=c(min(X, -3), max(X, 3)))
x <- seq(-3, 3, 0.03)
lines(x, dnorm(x))

# visualisation
sectionview(m, center = rep(0.5, 5))

#### Question 10 ####
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
