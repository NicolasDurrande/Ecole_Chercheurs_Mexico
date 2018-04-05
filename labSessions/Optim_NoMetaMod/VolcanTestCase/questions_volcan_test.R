#######################################################
# OPTIMIZATION WITHOUT METAMODEL - PRACTICAL SESSION
# WITH THE VOLCANO TEST CASE
#######################################################
#?????????????????????????????????????????????
# Version avec questions en francais sans accents.
# Qxxx: ceci est une question?
#?????????????????????????????????????????????


rm(list=ls()) #  cleaning up
.defaultpar <- par() # useful to restaure default sometimes par(.pardefault)

#######################################################
#   VOLCANO TEST CASE PRESENTATION
#######################################################
# cf. documentation in the folder, file volcan_test_case.pdf
#
source("utilities_volcan.R")
#
# The volcano is a calibration problem in 5 continuous variables:
# The parameters are normalized. The solution that we will seek through 
# optimization later is 
cat("normalized xstar =",norm_var(xstar))
cat("the number of variables (d in class, nbvar in the code) is ",nbvar)
#
# a random "volcano" with normalized model-measure distance nwls
a_normalized_volc <-runif(nbvar) 
nwls <- compute_wls(a_normalized_volc)
cat("nwls=",nwls,"\n")
# denormalized associated characteristics of the volcano
cat("the volcano is :\n")
print(unnorm_var(a_normalized_volc))
# details about the ground displacements associated to this volcano 
# could be plotted as contours by following example in the file plots_3d_full_grid.R
# xsol1 <- unnorm_var(a_normalized_volc)

#######  design of experiments #############################
library(DiceDesign)
nbinit <- 100 # number of points in the initial design of experiments
set.seed(42)

# do an optimal Latin Hypercube Sampling
Xnorm <- lhsDesign(n=nbinit,dimension=nbvar)$design
colnames(Xnorm) <- varnames
pairs(Xnorm) # look at it
# calculate weighted least squares at X of U projected on LOS (w.r.t. target Glb_ulos)
norm_wls <- compute_wls(Xnorm)
hist(norm_wls)
# for the function distribution, the classes of points that are
# good, pretty good, average, bad, are separated by the levels
zlevels<-c(-1.8,-1.2,0.3)

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
  res_rand <- random_search(compute_wls,param)
  rand_fhist[,i] <- res_rand$fhist
  # plot optimization results, there are at least 2 plots done
  # notice that, of course, there are variations in the results (it's random)
  colnames(res_rand$xhist)<-varnames
  #   in the pairs plots, the order of colors for good (low f)
  #   to bad (high f) points is: red, orange, lightblue, lightgreen
  resopt_plot(res_rand$fhist,res_rand$xhist,res_rand$x_best,res_rand$f_best,
              flevels=zlevels,log=F)
  cat("random_search run ",i," done")
  readline(prompt="Press [enter] to continue")
}
# note on the pairs plots that there seems to be a continuum 
# of red dots (good solutions) particularly in the (a,p) plane.
# Indeed, there is an infinity of solutions in this problem that 
# are associated to a relation a^3*p = astar^3*pstar (cf. mogi_3D.R)

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

#?????????????????????????????????????????????
# Q9
# brancher la fonction qui calcule la distance moindre carres 
# du volcan sur la fonction globale qui est optimisee dans le 
# reste du code, i.e., remplir
#  fun <<- 
#?????????????????????????????????????????????

#
########
# look at the effect of the initial point
init_glob_record_f()
x0 <- rep(0.75,nbvar) # an initial point
opts <- list(xtol_rel = 1e-14, maxeval = budget) # options for Nelder Mead
res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_nm$par,res_nm$value,flevels = zlevels,log = F)
# this seems to be a nice case of Nelder-Mead crashing on a variable bound
# and subsequently degenerating

# 2nd try
init_glob_record_f()
x0 <- rep(0.25,nbvar) # an initial point
opts <- list(xtol_rel = 1e-14, maxeval = budget) # options for Nelder Mead
res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_nm$par,res_nm$value,flevels = zlevels,log = F)
# which converges to the true optimum
cat("distance between NM final point and opt =",norm(x = (res_nm$par-norm_var(xstar)),type="2"))

# 3rd try: do an early stop AND ATTEMPT A RESTART
init_glob_record_f()
x0 <- rep(0.75,nbvar) # an initial point
# set stop option at a large simplex size
opts <- list(xtol_rel = 1e-2, maxeval = budget) 
res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
# restart from convergence point, remove used budget
x0 <- res_nm$par
opts <- list(xtol_rel = 1e-14, maxeval = (budget-length(glob_fhist))) # options for Nelder Mead
res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_nm$par,res_nm$value,flevels = zlevels,log = F)
# look at the pairs plot. Here the first part of the run degenerates on 3 bounds, 
# which is too much for a restart to save it.
cat("distance between NM final point and opt =",norm(x = (res_nm$par-norm_var(xstar)),type="2"))

########
# repeated Nelder-Mead runs from many random points to have
# a more complete understanding of how it behaves on this problem
nm_fhist <- matrix(,nrow = budget,ncol = no_test)
nm_xbest_hist <- matrix(,nrow=no_test,ncol=nbvar)
for(i in 1:no_test){
  init_glob_record_f()
  x0 <- runif(n = nbvar,min = 0,max = 1) # starting point
  opts <- list(xtol_rel = 1e-14, maxeval = budget) # some options for the optimizer
  res_nm <- neldermead(x0,ofwrapper, lower = rep(0,nbvar), upper = rep(1,nbvar),control = opts)
  if (res_nm$value < -2.2) {nm_xbest_hist[i,]<-res_nm$par}
  store_hist <<- FALSE # end history recording
  if (nrow(glob_fhist) >= budget) {
    nm_fhist[,i] <- glob_fhist[1:budget,1]}
  else {
  nm_fhist[,i] <- rbind(glob_fhist,matrix(NA,nrow = budget - nrow(glob_fhist),ncol = 1))  
  }
}
# interesting to see how the best points are spread
nm_xbest_hist<- nm_xbest_hist[which(!is.na(nm_xbest_hist[,1])),] # remove bad runs
colnames(nm_xbest_hist)<-varnames
pairs(nm_xbest_hist,xlim=c(0,1),ylim=c(0,1))
title("location of good NM convergence points",line = +3)

# compare results of Nelder-Mead with those of random_search
all_hist <- cbind(rand_fhist,nm_fhist)
comp_opt_plot(all_hist=all_hist,no_test=no_test,lab_opt=c("random search","Nelder-Mead"),pos="bottomright",log=F)
# random_search has more spread in performance than in the tension
# test case which is an illustration that the function f is not as regular.
# Nelder-Mead converges in about 50% of the cases, 
# and fails at a level worse than random search the rest of the time. 
# Relatively, random search is robust.
# As seen in the detailed runs above, it is likely that the variables 
# bounds create this false convergence.
# Nelder-Mead typically does not use all the budget.


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
x0 <- rep(0.5,nbvar) # an initial point
opts_defaults <- opts <- list(sigma=0.03,lambda=zlambda,mu=zmu,vectorized=F,
             maxit=ceiling(budget/zlambda),stopfitness=-Inf)
res_cma <- cma_es(par = x0, ofwrapper,lower=rep(0,nbvar), 
                  upper=rep(1,nbvar), control=opts)
colnames(glob_xhist)<-varnames
par(.defaultpar)
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,flevels = zlevels,log = F)
# works great. I had to decrease the step size wrt to the tension test
# case to get it to work. I think the way the bounds are handled in cmaes
# is bad and ruins convergence, and the bounds are hit too often if 
# sigma is large

# again with another starting point
x0 <- rep(0.75,nbvar)
init_glob_record_f()
res_cma <- cma_es(par = x0, ofwrapper,lower=rep(0,nbvar), 
                  upper=rep(1,nbvar), control=opts)
colnames(glob_xhist)<-varnames
par(.defaultpar)
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,flevels = zlevels,log = F)
# bad convergence because the bounds are hit ... 

# start with a smaller step size
opts$sigma <- 0.01
x0 <- rep(0.75,nbvar)
init_glob_record_f()
res_cma <- cma_es(par = x0, ofwrapper,lower=rep(0,nbvar), 
                  upper=rep(1,nbvar), control=opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,flevels = zlevels,log = F)
# the smaller step fixes the bound issue ... 

# finally try with a large population
par(.defaultpar)
opts$sigma <- 0.01 # I keep a really small step size, otherwise does not work
opts$lambda <- 300
opts$mu <- 150
opts$maxit=ceiling(budget/opts$lambda)
x0 <- rep(0.75,nbvar)
init_glob_record_f()
res_cma <- cma_es(par = x0, ofwrapper,lower=rep(0,nbvar), 
                  upper=rep(1,nbvar), control=opts)
colnames(glob_xhist)<-varnames
resopt_plot(glob_fhist,glob_xhist,res_cma$par,res_cma$value,flevels = zlevels,log = F)
# larger population slows down the search.
# Notice that the valley of optimal points is visited in part (red 
# curve in the (a,p) plane)

####
# do repeated tests for comparison with the other optimizers
cma_fhist <- matrix(,nrow = budget,ncol = no_test)
cma_xbest_hist <- matrix(,nrow=no_test,ncol=nbvar)
opts_defaults$sigma <- 0.01 # small is beautiful with this implementation...
for (i in 1:no_test) {
  init_glob_record_f()
  res_cma <- cma_es(runif(n=nbvar), ofwrapper,lower=rep(0,nbvar), 
                    upper=rep(1,nbvar), control=opts_defaults)
  if (res_cma$value < -2.2) {cma_xbest_hist[i,]<-res_cma$par}
  if (nrow(glob_fhist) >= budget) {
  cma_fhist[,i] <- glob_fhist[1:budget,1]}
  else {
    cma_fhist[,i] <- rbind(glob_fhist,matrix(NA,nrow = budget - nrow(glob_fhist),ncol = 1))  
  }
}
# and plot
# compare results of Nelder-Mead with those of random_search
all_hist <- cbind(rand_fhist,nm_fhist,cma_fhist)
comp_opt_plot(all_hist=all_hist,no_test=no_test,lab_opt=c("random search","Nelder-Mead","CMA-ES"),pos="nolegend",log=F)
# interesting to see how the best points are spread
cma_xbest_hist<- cma_xbest_hist[which(!is.na(cma_xbest_hist[,1])),] # remove bad runs
colnames(cma_xbest_hist)<-varnames
pairs(cma_xbest_hist,xlim=c(0,1),ylim=c(0,1))
title("location of good CMA convergence points",line = +3)

#?????????????????????????????????????????????
# Q10
# Analyser les comparaisons entre optimiseurs sur le cas test volcan
#?????????????????????????????????????????????




# Like with the tension test case, Nelder Mead fast and unreliable, 
# cma-es a bit slower but more accurate. Yet it will often get stuck
# at the bounds. There might be an issue with bounds handling 
# in this implementation of cma-es...
# Notice that the set of convergence points cma_xbest_hist clearly 
# shows the indetermination in a and p (the other variables are well fixed)

