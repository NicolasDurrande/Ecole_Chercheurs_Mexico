library(DiceKriging)
library(DiceView)

#### Question 1 ####
load('toy_data.Rdata')
plot(X, Y)

#### Question 2 ####
m <- km(y~1, design=data.frame(x=X), response=data.frame(y=Y), 
        covtype="gauss", coef.trend = 0, coef.var=1, coef.cov=.2)

x <- matrix(seq(-0.5, 1.5, 0.01))
Z <- simulate(m, 10, newdata=data.frame(x=x), cond=FALSE, nugget.sim=1e-8)
Z <- t(Z)

matplot(x, Z, type='l', col=1)

#### Question 3 ####
Zc <- simulate(m, 10, newdata=data.frame(x=x), cond=TRUE, nugget.sim=1e-8)
Zc <- t(Zc)
matplot(x, Zc, type='l', col=1)
points(X, Y)

sectionview(m, xlim=c(-0.2, 1.2), title='GP regression')

leaveOneOut.km(m, type='UK')

