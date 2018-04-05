#######################################################
# OPTIMIZATION WITHOUT METAMODEL - PRACTICAL SESSION
# WITH THE TENSION TEST CASE
#
#?????????????????????????????????????????????
# Version avec questions en francais sans accents.
# Qxxx: ceci est une question?
#?????????????????????????????????????????????
#
#######################################################

rm(list=ls()) #  cleaning up
.defaultpar<-par()
#######################################################
#   TENSION TEST CASE PRESENTATION
#######################################################
# load the test case. Provides the function to minimize in 4 variables
# j2_tension_norm(nt = x)
# where nt is a normalized (between 0 and 1) vector of variables
#
source('./tension_test.R')
#
# The problem "tension test" is a calibration problem:
# find the 4 parameters of a metal law so that the results of a tension 
# test match an "experimental" tension test in a least squares sense. 
# cf. https://en.wikipedia.org/wiki/Stress%E2%80%93strain_curve
# cf. https://en.wikipedia.org/wiki/Work_hardening
# The 4 parameters are Young's modulus, the elasticity limit, 
# the hardening speed and the hardening deformation. 
# The parameters are normalized between 0 and 1. 
# The solution that we will seek through optimization later is 
cat("normalized xstar =",norm_var(x=theta_star))
# note that the true variables are normalized, which is important here
# because the variables have orders of magnitude differences: typically
# Young's moduli are in the 100000's MPa and the hardening deformation
# in the 0.001. The functions norm_var() and unnorm_var() implement
# normalization of the variables between 0 and 1 and unnormalization.
cat("the number of variables (d in class, nbvar in the code) is ",nbvar)
#
# Let us first generate a few candidate solutions and plot the associated 
# tension curve, the target tension curve, and the associated J2 distance
#
ntry <- 4 # nb of random candidate solutions (materials here)
for (i in 1:ntry) {
  ntheta <- runif(1,min=0.2,max=1.8)*norm_var(theta_star) # in calibration, x often called theta (here normalized)
  # Note: if ntheta <- runif(nbvar,min=0,max=1) the curves are so far from target that it is hard to plot
  sig <- tension_test(eps_tot = eps_tot,theta = unnorm_var(ntheta))
  j2 <- j2_tension_norm(nt = ntheta)
  # plot the target curve, all defined in the file "tension_test.R"
  plot(eps_tot,sig_star,type="o",col="black",ylim = c(0,600),xlab = "total strain",ylab = "stress")
  # superimpose the random material
  lines(eps_tot,sig,type = "o",col="red")
  legend(x = "bottomright",legend = c("target","candidate"),col=c("black","red"),lty=c(1,1))
  t1 <- paste(format(ntheta,scientific = T),collapse =" ")
  t1 <- paste("theta = ",t1,"\nj2 = ",j2)
  title(t1,cex.main=0.8)
}
#
# I like this test case because it is rooted in real world problems:
# it has the flavor of true calibration problems. J2(theta) is continuous
# but not continuously differentiable because of a thresholding effect 
# between the elastic and plastic domains of deformation, and J2 is not 
# very well conditioned (I guess).
# Look at the objective function distribution with a design of experiments 
library(DiceDesign)
nbinit <- 100 # number of points in the initial design of experiments
set.seed(42) # for repeatability
# do an optimal Latin Hypercube Sampling
Tnorm <- lhsDesign(n=nbinit,dimension=nbvar)$design
colnames(Tnorm) <- varnames
pairs(Tnorm) # look at it
# calculate weighted least squares at the thetas
j2s <- apply(X = Tnorm,MARGIN = 1,FUN = j2_tension_norm)
# look at the function distribution: 
hist(j2s,breaks = 10) 
# j2_tension_norm has large values, is clearly not gaussian so modeling
# it with a Gaussian process would probably require further normalization 
# on J2 this time.

#######################################################
#  REPEATED IDENTIFICATIONS WITH A RANDOM OPTIMIZER
#######################################################
#  Always a good idea to compare an optimizer with a random search:
#  a random sampling is insensitive to every mathematical features of
#  the optimized function but it is really slow (pays the full price
#  of the curce of dimensionality).
#  Random searches have to be repeated to get a true feeling 
# (away from lucky/unlucky results)
no_test <- 10 # number of repetitions of the optimization run
budget <- 4000 # budget of each optimization (nb calls to the objective function)
# load a provided random optimizer, ...
#   really just wrapping around runif(nbvar,min=0,max=1)
source('../Optimizers/RSalgorithm.R')
# load a function to plot results of optimization runs
source('../Optimizers/resopt_plot.R')
# initialize a matrix to save search results for later comparisons.
# typical syntax is "name_of_optim"_fhist for storing history of function
rand_fhist <- matrix(,nrow = budget,ncol = no_test)
# parameters to random_search, hopefully self explicit
param <- list(LB=rep(0,nbvar),UB = rep(1,nbvar),budget = budget,dim=nbvar)
for(i in 1:no_test){
  res_rand <- random_search(j2_tension_norm,param)
  rand_fhist[,i] <- res_rand$fhist
  # plot optimization results, there are at least 2 plots done
  # notice that, of course, there are variations in the results (it's random)
  colnames(res_rand$xhist)<-varnames
  #   in the pairs plots, the order of colors for good (low f)
  #   to bad (high f) points is: red, orange, lightblue, lightgreen
  resopt_plot(res_rand$fhist,res_rand$xhist,res_rand$x_best,res_rand$f_best,c(100,1000,10000))
  cat("random_search run ",i," done")
  readline(prompt="Press [enter] to continue")
}
#?????????????????????????????????????????????
# Q1: en regardant les pairs plots, 
#?????????????????????????????????????????????

########################################################
# OPTIMIZE WITH NELDER MEAD
########################################################
# implementation from the nlopt library
# Documentation through ?neldermead and ?nl.opts
library("nloptr")
# Global variables (nbcalls, glob_fhist, glob_hist) 
# are used to store the optimization history by wrapping
# the objective function, cf. of_wrapper.R . This is a nice generic 
# way to monitor what packages do when you don't have access to their 
# sources.
# Initialize these global variable by just calling init_glob_record_f()
source("../Optimizers/of_wrapper.R")
# define the function to optimize for the wrapper
fun <<- j2_tension_norm
#
########
# look at the effect of the initial point
init_glob_record_f()
x0 <- c(0.75,0.75,0.75,0.75) # an initial point
opts <- list(xtol_rel = 1e-14, maxeval = budget) # options for Nelder Mead
res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_nm$par,res_nm$value,c(100,1000,10000))
# this seems to be a nice case of Nelder-Mead crashing on a variable bound
# and subsequently degenerating

#?????????????????????????????????????????????
# Q2: faire une autre optimisation par Nelder-Mead 
# en partant de (0.25,0.25,0.25,0.25) avec les memes options.
# L'optimiseur converge-t-il vers la solution?
#?????????????????????????????????????????????
# 2nd try

# ???? (...)

colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_nm$par,res_nm$value,c(100,1000,10000))
cat("distance between NM final point and opt =",norm(x = (res_nm$par-norm_var(x = theta_star)),type="2"))

#?????????????????????????????????????????????
# Q3: faire une double optimisation par Nelder-Mead 
# en partant de (0.75,0.75,0.75,0.75) avec
# une premiere optimisation qui s'arrete tot 
# (utiliser xtol_rel = 1e-3 ou plus grand)
# Puis redemarrer depuis le point de convergence 
# avec un autre Nelder-Mead ou xtol_rel est remis 
# a une valeur faible (1.e-14) et le budget ajuste 
# a ce qui reste 
# opts$maxeval = (budget-length(glob_fhist))
# L'optimiseur converge-t-il vers la solution?
#?????????????????????????????????????????????
# 3rd try: do an early stop AND ATTEMPT A RESTART
init_glob_record_f()
x0 <- c(0.75,0.75,0.75,0.75) # an initial point

# ???? (...)

colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_nm$par,res_nm$value,c(100,1000,10000))
cat("distance between NM final point and opt =",norm(x = (res_nm$par-norm_var(x = theta_star)),type="2"))

########
# repeated Nelder-Mead runs from many random points to have
# a more complete understanding of how it behaves on this problem
nm_fhist <- matrix(,nrow = budget,ncol = no_test)
for(i in 1:no_test){
  init_glob_record_f()
  x0 <- runif(n = nbvar,min = 0,max = 1) # starting point
  opts <- list(xtol_rel = 1e-14, maxeval = budget) # some options for the optimizer
  res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
  # plot(1:nbcalls,log10(glob_fhist))
  store_hist <<- FALSE # end history recording
  if (nrow(glob_fhist) >= budget) {
    nm_fhist[,i] <- glob_fhist[1:budget,1]}
  else {
  nm_fhist[,i] <- rbind(glob_fhist,matrix(NA,nrow = budget - nrow(glob_fhist),ncol = 1))  
  }
}

# compare results of Nelder-Mead with those of random_search
all_hist <- cbind(rand_fhist,nm_fhist)
par(.defaultpar)
comp_opt_plot(all_hist=all_hist,no_test=no_test,lab_opt=c("random search","Nelder-Mead"),pos="bottomright")
# Note that the y-axis has a log10 scale.

#?????????????????????????????????????????????
# Q4: que peut-on dire de la comparaison entre 
# Nelder-Mead et random search?
#?????????????????????????????????????????????

########################################################
# OPTIMIZE WITH CMA-ES
########################################################
# implementation from the cmaes library
# Documentation through ?cma_es
library("cmaes")
# I'm not reloading the of_wrapper and redefining fun, it was done with Nelder-Mead

####
# do a few single runs changing the initial point, initial step size
# and population sizes
init_glob_record_f()
# initialize parameters of CMA-ES
zlambda <- 4+floor(3*log(nbvar)) # number of samples, aka population size
zmu <- floor(zlambda/2) # number of points kept in update
x0 <- c(0.75,0.75,0.75,0.75) # an initial point
opts_defaults <- opts <- list(sigma=0.3,lambda=zlambda,mu=zmu,vectorized=F,
             maxit=ceiling(budget/zlambda),stopfitness=-Inf)
res_cma <- cma_es(par = x0, ofwrapper,lower=rep(0,nbvar), 
                  upper=rep(1,nbvar), control=opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,c(100,1000,10000))
# usually works well, slower than Nelder-Mead


#?????????????????????????????????????????????
# Q5: faire une autre optimisation par CMA-ES
# en partant de (0.25,0.25,0.25,0.25) avec les memes options.
# L'optimiseur converge-t-il vers la solution?
#?????????????????????????????????????????????
# 2nd try

# ???? (...)

colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,c(100,1000,10000))

#?????????????????????????????????????????????
# Q6: faire une autre optimisation par CMA-ES
# en partant de (0.75,0.75,0.75,0.75) avec une 
# taille de pas petite, de l'ordre de 0.01
# Que remarquez vous?
#?????????????????????????????????????????????
init_glob_record_f()

#???? (...)

colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,c(100,1000,10000))
# the convergence is longer than the other one started from the same point

#?????????????????????????????????????????????
# Q7: faire une autre optimisation par CMA-ES
# en partant de (0.75,0.75,0.75,0.75) avec une 
# population (lambda) grande. Penser a ajuster 
# mu et maxit, remettre une taille de pas classique.
# Que remarquez vous?
#?????????????????????????????????????????????
# finally try with a large population
init_glob_record_f()

#???? (...)

colnames(glob_xhist)<-varnames
par(.defaultpar)
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,c(100,1000,10000))

####
# do repeated tests for comparison with the other optimizers
cma_fhist <- matrix(,nrow = budget,ncol = no_test)
for (i in 1:no_test) {
  init_glob_record_f()
  res_cma <- cma_es(runif(n=nbvar), ofwrapper,lower=rep(0,nbvar), upper=rep(1,nbvar), control=opts_defaults)
  if (nrow(glob_fhist) >= budget) {
  cma_fhist[,i] <- glob_fhist[1:budget,1]}
  else {
    cma_fhist[,i] <- rbind(glob_fhist,matrix(NA,nrow = budget - nrow(glob_fhist),ncol = 1))  
  }
}
# and plot
# compare results of Nelder-Mead with those of random_search
all_hist <- cbind(rand_fhist,nm_fhist,cma_fhist)
comp_opt_plot(all_hist=all_hist,no_test=no_test,lab_opt=c("random search","Nelder-Mead","CMA-ES"),pos="nolegend")

#?????????????????????????????????????????????
# Q8: Comparer les comportements des trois optimiseurs
# (random, Nelder-Mead et CMA-ES). 
#?????????????????????????????????????????????

