library(DiceKriging)
library(DiceView)

#### Question 1 ####

#### Question 2 ####
m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y),
        coef.trend = 0, nugget=1e-8,
        covtype="gauss", coef.var=4, coef.cov=.5)

x <- matrix(seq(-0.2, 1.2, 0.01))
Z <- simulate(m, 10, newdata=data.frame(x=x), cond=FALSE)
Z <- t(Z)

matplot(x, Z, type='l', col=1)

#### Question 3 ####

#### Question 4 ####

#### Question 5 ####

res_std <- ... # à vous de jouer !

hist(res_std, freq=FALSE, xlim=c(min(X, -3), max(X, 3)))
x <- seq(-3, 3, 0.03)
lines(x, dnorm(x))

#### Question 6 ####

#### Question 7 ####

#### Question 8 ####

#### Question 9 ####

# test de la moyenne
res <-leaveOneOut.km(m, type='SK')
Q2 <- 1 - sum((Y - res$mean)^2) / sum((Y - mean(Y))^2)

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