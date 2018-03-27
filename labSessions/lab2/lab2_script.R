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
