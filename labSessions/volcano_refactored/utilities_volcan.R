###########################
# Utility functions for
# Inversion of a punctual displacements source from 3D data
# The data used are under-sampled with the quadtree method (irregular grid)
#
# Rodolphe Le Riche, Nicolas Durrande, Valerie Cayol, Victor Picheny,  Sebastien Roux
###########################

####### load utilities ##########################
mydist <- function(x,y,theta=NULL){
  if(ncol(x) != ncol(y)) stop("x and y must have the same number of columns") 
  if(is.null(theta)) theta <- rep(1,ncol(x))
  dist2 <- matrix(0,dim(x)[1],dim(y)[1])
  for(i in 1:dim(x)[2]){
    dist2 <- dist2 + (outer(x[,i],y[,i],"-")/theta[i])^2
  }
  return(sqrt(dist2))
}

kExp <- function(x,y=NULL,covparam=NULL){
  if (is.null(covparam)) covparam <- c(1,rep(.2,ncol(x)))
  if (is.null(y)) y <- x
  return(abs(covparam[1])*exp(-mydist(x,y,covparam[-1])))
}

###### input for variables identification ###########
# xs: X location of source in m UTM coordinates
# ys: Y location of source in m UTM
# zs: Elevation of source with respect to sea level in m
# a: source radius in m
# p: Source overpressure in MPa
varnames <- c("xs","ys","zs","a","p")
nbvar <- 5

# optimum and bounds on variables
xstar <- c(367000, 7650300, 0, 500, 20)
xmax  <- c(368000, 7651000, 1000, 1000, 500)
xmin  <- c(364000, 7649000, -3000, 100, -500)
names(xstar) <- names(xmax) <- names(xmin) <- varnames
Glb_var <- list(nbvar=nbvar, xmax=xmax, xmin=xmin) # always useful stuff

#####################################################################################
####### load data ###################################################################
#####################################################################################
data <- read.csv(file="data_nonoise.csv")
Glb_xi <- as.matrix(data[1])
Glb_yi <- as.matrix(data[2])
Glb_zi <- as.matrix(data[3])
Glb_ulos <- as.matrix(data[4])

# calculate data Covariance matrix, store it in a Global variable
# covariance from exponential kernel, var = 5e-4m2, cor_length = 850 m
# and invert it
Glb_CXinv <- solve(kExp(x=cbind(Glb_xi, Glb_yi), covparam=c(5e-4,850,850))) # calculated once for all, used in wls_ulos

Glb_data <- list(xi=Glb_xi, yi=Glb_yi, zi=Glb_zi, ulos=Glb_ulos, CXinv=Glb_CXinv)

#####################################################################################
mogi_3D <- function(G, nu, xs, ys, zs, a, p, xi, yi, zi){
  # MOGI(G,nu,xs,ys,zs,a,p,xi,yi,zi) compute surface displacements and tilts created by
  # a point source located beneath a topography. To account for topography, a
  # first order solution in which the actual source to ground surface point
  # is taken into account
  # 
  # Parameters are 
  # G = shear modulus in MPa, G = E/2(1+nu)
  # nu = Poisson's ratio
  # xs, ys, zs = source position (z axis is positive upward),
  # a = source radius
  # p = source overpressure in MPa, 
  # xi, yi, zi = location of ground surface points
  #
  # V. Cayol, LMV, sept 2017
  # (translated into R by R. Le Riche)
  
  DV <- pi*a^3*p/G
  C <- (1-nu)*DV/pi
  r <- sqrt((xi-xs)^2+(yi-ys)^2)
  f <- r^2+(zi-zs)^2
  uzi <- C*(zi-zs)/(f^(3/2))
  ur <- C*r/(f^(3/2))
  theta <- atan2(yi-ys,xi-xs)
  uxi <- ur*cos(theta)
  uyi <- ur*sin(theta)
  return(list(x=uxi,y=uyi,z=uzi))
}
#####################################################################################
wls_ulos <- function(xx, Glb_data){
  # Weighted Least Squares distance function for ulos vectors
  # xx must be a vector of size 5 with names (xs, ys, zs, a, p)
  # Glb_data is a list that contains {xi, yi, zi, ulos, CXinv}
  
  G = 2000 # Shear modulus in MPa
  nu = 0.25 # Poisson's ratio
  nlos = c(-0.664,-0.168,0.728) # vector of direction of line of sight (satellite)

  # Compute surface displacements
  U <- mogi_3D(G=G, nu=nu, xs=xx["xs"], ys=xx["ys"], zs=xx["zs"], a=xx["a"], p=xx["p"],
               xi=Glb_data$xi, yi=Glb_data$yi, zi=Glb_data$zi)
  
  # project along los
  ulos <- nlos[1]*U$x+nlos[2]*U$y+nlos[3]*U$z
  
  # calculate weighted least squares (WLS)
  return(t((ulos-Glb_data$ulos))%*%Glb_data$CXinv%*%(ulos-Glb_data$ulos))
}

#####################################################################################
#######  useful functions ###########################################################
#####################################################################################
# Scale from x in [0 1] to xx in [xmin xmax]
unnorm_x <- function(x){
  if (is.null(dim(x))) x <- matrix(data = x, nrow=1) # numeric vector
  nbrep <- nrow(x)
  xx <- matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar) + 
    x * matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)
  colnames(xx) <- varnames
  return(xx)
}
#####################################################################################
# Scale from xx in [xmin xmax] to x in [0 1] 
norm_x <- function(xx){
  if (is.null(dim(xx))) xx <- matrix(data = xx, nrow=1)
  nbrep <- nrow(xx)
  x <- (xx - matrix(rep(xmin,times=nbrep),byrow = T,ncol=nbvar)) / 
    matrix(rep((xmax-xmin),times=nbrep),byrow = T,ncol=nbvar)  
  colnames(x) <- varnames
  return(x)
}
#####################################################################################
# Normalized weighted least square as a function of a vector in [0,1]^5
compute_wls <- function(x) {
  # x can be a vector or a matrix with values in [0,1]
  if (is.null(dim(x))) {
    x <- matrix(data = x, nrow=1)
  }  else if (ncol(x)!=nbvar) {
    x <- t(x)
  }
  xx <- unnorm_x(x=x)
  y <- apply(xx, 1, wls_ulos)
  # normalize the output so that it is centered with a unit std dev
  # because wls ranges from 0 to 10^9, do a scaling in log(1+wls)
  return((log(1+y) - 8.54)/3.2)
}
#####################################################################################
# Norm of the displacement vector at a location=(x,y,z) given given parameters xx
displacement <- function(xx, location){
  # xx must be a vector of size 5 with names (xs, ys, zs, a, p)
  # location must be a vector of size 3 (xi, yi, zi)
  
  G <- 2000 # Shear modulus in MPa
  nu <- 0.25 # Poisson's ratio
  nlos <- c(-0.664,-0.168,0.728) # vector of direction of line of sight (satellite)
  
  # Compute surface displacements at location 
  U <- mogi_3D(G=G, nu=nu, xs=xx["xs"], ys=xx["ys"], zs=xx["zs"], a=xx["a"], p=xx["p"],
               xi=location[1], yi=location[2], zi=location[3])
  
  # norm of displacement (U$x*U$x + U$y*U$y + U$z*U$z)
  u_norm <- sqrt(sum(unlist(U)^2))
  names(u_norm) <- "distance"
  return(u_norm)  
}
#####################################################################################
# Norm of the displacement vector at a location=(x,y,z) as a function of a vector in [0,1]^5
compute_displacement <- function(x, location) {
  # x can be a vector or a matrix with values in [0,1]
  # location must be a vector of size 3 (xi, yi, zi)

  if (is.null(dim(x))) {
    x <- matrix(data = x, nrow=1)
  }  else if (ncol(x)!=nbvar) {
    x <- t(x)
  }
  xx <- unnorm_x(x=x)
  y <- apply(xx, 1, displacement, location)
  return(y)
}