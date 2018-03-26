
f <- function(x){
  2 + 3*x^2 - sin(6*x)
}

X <- runif(10)
Y <- f(X)

plot(X, Y)

save(X, Y, file='toy_data.Rdata')

load('toy_data.Rdata')
