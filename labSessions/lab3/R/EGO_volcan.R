######### EXAMPLE OF USE OF THE VOLCANO TEST CASE WITH KRIGING MODEL ###########
######### BUILDING AND EGO OPTIMIZATION                              ###########
source("utilities_volcan.R")

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
plot.ecdf(norm_wls) # look at it

# plot learning data
par(mfrow=c(2,3))
for (i in 1:nbvar)  plot(Xnorm[,i],norm_wls,xlab=colnames(Xnorm)[i])
# empty plot + text
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("LEARNING\n","SET"), cex = 1.6, col = "black")

# Build initial kriging model
library(DiceKriging)
mod <- km(~1, design=Xnorm, response=norm_wls, control=list(trace=FALSE), optim.method = "gen",
          lower=rep(.1,5), upper=rep(1,5))

# Run EGO
library(DiceOptim)
res <- EGO.nsteps(model=mod, fun=compute_wls, nsteps=20, lower=rep(0, 5), upper=rep(1,5), control=list(print.level=0))
xbest <- res$lastmodel@X[which.min(res$lastmodel@y),]

# last model
resLOO <- leaveOneOut.km(res$lastmodel, type="UK", trend.reestim=FALSE)
Q2 <- 1 - sum((resLOO$mean - res$lastmodel@y)^2) / sum( (res$lastmodel@y - mean(res$lastmodel@y))^2)

# Draw results
par(mfrow=c(1,2))
hist(norm_wls)
hist(res$value)
plot(res$lastmodel)

library(DiceView)
sectionview(model=res$lastmodel, center=xbest)

cols=c(rep("blue", 100), rep("red", 50))
cols[which(res$lastmodel@y < -1.6)] <- "green"
X2 <- res$lastmodel@X
pairs(X2, col=cols)
