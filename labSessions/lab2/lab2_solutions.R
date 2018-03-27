library(DiceKriging)
library(DiceView)

#### Question 1 ####
load('toy_data.Rdata')
plot(X, Y)

#### Question 2 ####
m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y),
        coef.trend = 0, nugget=1e-8,
        covtype="gauss", coef.var=4, coef.cov=.5)

x <- matrix(seq(-0.2, 1.2, 0.01))
Z <- simulate(m, 10, newdata=data.frame(x=x), cond=FALSE)
Z <- t(Z)

matplot(x, Z, type='l', col=1)

#### Question 3 ####
Zc <- simulate(m, 10, newdata=data.frame(x=x), cond=TRUE)
Zc <- t(Zc)
matplot(x, Zc, type='l', col=1)
points(X, Y)

sectionview(m, xlim=c(-0.2, 1.2), title='GP regression')

#### Question 4 ####

m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", coef.trend = 0, nugget=1e-8)

print(m)
sectionview(m, xlim=c(-0.2, 1.2), title='GP regression')

#### Question 5 ####
res <-leaveOneOut.km(m, type='SK')
res_std <- (Y-res$mean) / res$sd

hist(res_std, freq=FALSE, xlim=c(min(X, -3), max(X, 3)))
x <- seq(-3, 3, 0.03)
lines(x, dnorm(x))
plot(m)

#### Question 6 ####
predict.km(m, newdata=data.frame(x=10), type='SK', light.return = TRUE)

sectionview(m, xlim=c(-0.2, 4), ylim=c(-18, 22), title='GP regression')

# modèle de krigeage ordinaire (estimation tendance constante)
m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", nugget=1e-8)
sectionview(m, xlim=c(-0.2, 4), ylim=c(-18, 22), title='GP regression')

# modèle de krigeage universel (estimation tendance lineaire)
m <- km(y~x, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", nugget=1e-8)
sectionview(m, xlim=c(-0.2, 4), ylim=c(-18, 22), title='GP regression')

#### Question 7 ####
load('XY_volcano.Rdata')

# test des données chargées
par(mfrow=c(1,5))
for(i in 1:5){
  plot(X[,i], Y)
}
par(mfrow=c(1,1))

#### Question 8 ####

m <- km(y~., design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="matern5_2", nugget=1e-8)

print(m)

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
pdf("sobolGP.pdf")
plot(res)
dev.off()
