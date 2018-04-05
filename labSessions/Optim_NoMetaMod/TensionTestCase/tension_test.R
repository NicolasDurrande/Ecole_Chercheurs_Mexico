#####################################################
# the hidden files for the ...
# least squares identification of a (house) model 
# that resembles a tension curve for metals
# (Rodolphe Le Riche, 2018)
#####################################################

#### inputs 
  # model parameters
  n2i <- list(E=1, R=2, c=3, d=4) # name to index for variables
  varnames <- c("E","R","c","d")
  nbvar <- 4
  # optimum
  theta_star <- NA
  theta_star[n2i$E] <- 100000 # Young's modulus (e.g., in MPa)
  theta_star[n2i$R] <- 200 # elasticity limit (e.g., in MPa)
  theta_star[n2i$c] <- 100 # hardening speed 1
  theta_star[n2i$d] <- 0.0015 # hardening "deformation" (*E to get stress) 
  # lower and upper bounds for identification problems
  theta_min <- NA
  theta_min[n2i$E] <- 40000 
  theta_min[n2i$R] <- 10 
  theta_min[n2i$c] <- 10 
  theta_min[n2i$d] <- 0.0001 
  theta_max <- NA
  theta_max[n2i$E] <- 300000 
  theta_max[n2i$R] <- 500
  theta_max[n2i$c] <- 500
  theta_max[n2i$d] <- 0.01 
  # total deformation is controlled
  eps_min <- 1.e-7
  eps_max <- 1.e-1
  eps_tot <- seq(from = eps_min, to = eps_max, length.out = 100)
#### end inputs
# a simplistic plasticity model, 1D in tension
tension_test <- function(eps_tot,theta){
  RoE <- theta[n2i$R]/theta[n2i$E]
  eps_p <- -RoE + eps_tot - theta[n2i$d]*(1-exp(-theta[n2i$c]*(eps_tot-RoE)))
  # could be made more complex by adding -d*(1-exp(-c*(eps_tot-R/E))) terms with other d and c
  eps_p[eps_tot<RoE] <- 0
  sig <- theta[n2i$E]*(eps_tot-eps_p)
  return(sig)
}

# target or experimental result
sig_star <- tension_test(eps_tot = eps_tot,theta = theta_star)
# plot(eps_tot,sig_star)

# (un)normalization utility functions
norm_var <- function(x){
  nt <- (x - theta_min)/(theta_max-theta_min)
  return(nt)
}

unnorm_var <- function(xn){
  tt <- theta_min+xn*(theta_max-theta_min)
  return(tt)
}

# weighted least square distance between results
j2_tension <- function(theta){
  nexp <- length(eps_tot)
  sig <- tension_test(eps_tot = eps_tot,theta = theta)
  W <- diag(nexp)
  j2 <- 1/(2*nexp)*(sig-sig_star)%*%W%*%(sig-sig_star)
  return(j2)
}

# j2 with normalized input
j2_tension_norm <- function(nt){
  return(j2_tension(unnorm_var(nt)))
}

