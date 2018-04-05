#################################################
# utilities to plot results of optimization runs
#
# (for now I plot log10(f), more readable for convergence, 
# but will fail with negative functions)
#################################################

# plot results of an optimization run
#   in the pairs plots, the order of colors for good (low f)
#   to bad (high f) points is: red, orange, lightblue, lightgreen
resopt_plot <- function(fh,xh,xb,fb,flevels,log=T){
  if (!length(flevels)==3) stop("3 function levels must be given")
  nlevels <- 3
  #### usual f vs calls plot
  if (log)  {plot(1:length(fh),log10(cummin(fh)),type="l",xlab = "nb. of calls",ylab = "log10(f)")}
  else {plot(1:length(fh),cummin(fh),type="l",xlab = "nb. of calls",ylab = "f")}

  chbest <- paste(format(xb,scientific = T),collapse = " ")
  chbest <- paste("best x = ",chbest,"\n f =",fb)
  title(chbest,cex.main=0.8)
  #### pairs plot
  zi <- list()
  ni <- rep(0,(nlevels+1))
  zi[[1]]<-which(fh<=flevels[1])
  ni[1]<-length(zi[[1]])
  for (i in 2:nlevels){
    zi[[i]]<-which((flevels[i-1]<fh)&(fh<=flevels[i]))
    ni[i]<-length(zi[[i]])
  }
  zi[[nlevels+1]]<-which(fh>flevels[nlevels])
  ni[nlevels+1]<-length(zi[[nlevels+1]])
  cols<-rep("lightgreen",length(fh))
  if (ni[3]>0) {cols[(ni[4]+1):(ni[4]+ni[3])]<-"lightblue"}
  if (ni[2]>0) {cols[(ni[4]+ni[3]+1):(ni[4]+ni[3]+ni[2])] <-"orange"}
  if (ni[1]>0) {cols[(ni[4]+ni[3]+ni[2]+1):(length(fh))]<-"red"}
  zcex<-rep(0.5,length(fh))
  pairs(rbind(xh[zi[[4]],],xh[zi[[3]],],xh[zi[[2]],],xh[zi[[1]],]), col=cols,pch=20,cex=zcex)
  t1 <- paste(format(flevels,scientific = T),collapse =" ")
  t1 <- paste("X for f levels ",t1)
  title(t1,cex.main=0.8,line = +3)
}

# get the 3 quartiles
threequartiles <- function(x){
  return(quantile(x,probs = c(0.25,0.5,0.75),na.rm = T,names=F))
}

# plot for comparing several optimizers
comp_opt_plot <- function(all_hist,no_test,lab_opt,pos="bottomleft",log=T) {
  nbopt <- length(lab_opt)
  nol <- "nolegend"
  zcol <- c()
  for (i in 1:nbopt){ zcol <- c(zcol,rep(i,no_test))}
  if (log)  {matplot(log10(all_hist),type="l",lty=1,col=zcol, xlab = "nb. calls", ylab="log10(f)")}
  else {matplot(all_hist,type="l",lty=1,col=zcol, xlab = "nb. calls", ylab="f")}
  zcolleg<-1:nbopt
  if (pos!=nol) {legend(x = pos,legend = lab_opt,col=zcolleg,lty=rep(1,nbopt))}
  title("iterates")
  # same plot with cumulated min
  if (log) {matplot(log10(apply(all_hist,2,cummin)),type="l",lty=1,col=zcol, xlab = "nb. calls", ylab="log10(f)")}
  else {matplot(apply(all_hist,2,cummin),type="l",lty=1,col=zcol, xlab = "nb. calls", ylab="f")}
  if (pos!=nol) {legend(x = pos,legend = lab_opt,col=zcolleg,lty=rep(1,nbopt))}
  title("bests")
  # plot the 3 quartiles
  all_quart <- matrix(,nrow = nrow(all_hist),ncol = 3*nbopt)
  for (i in 1:nbopt){
    if (log) {all_quart[,(1+(i-1)*3):(i*3)] <- 
      t(apply(log10(apply(all_hist[,(1+(i-1)*no_test):(i*no_test)],2,cummin)),1,threequartiles))}
    else {all_quart[,(1+(i-1)*3):(i*3)] <- 
      t(apply(apply(all_hist[,(1+(i-1)*no_test):(i*no_test)],2,cummin),1,threequartiles))}
  }
  zcol <- c()
  for (i in 1:nbopt){ zcol <- c(zcol,rep(i,3))}
  if (log) {matplot(all_quart,type="l",lty=1,col=zcol, xlab = "nb. calls", ylab="log10(f)")}
  else {matplot(all_quart,type="l",lty=1,col=zcol, xlab = "nb. calls", ylab="f")}
  if (pos!=nol) {legend(x = pos,legend = lab_opt,col=zcolleg,lty=rep(1,nbopt))}
  title("quartiles")
}