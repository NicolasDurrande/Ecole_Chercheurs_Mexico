# wrapper to the objective functions
# used to store history in global variables
ofwrapper <- function(xx){
  xx <- t(as.matrix(xx))
  of <- apply(xx, 1, fun)
  if (store_hist){
    nbcalls <<- nbcalls + 1
    if (nbcalls==1){glob_xhist[nbcalls,]<<- xx
                    glob_fhist[nbcalls]<<-of}
    else {glob_xhist<<- rbind(glob_xhist,xx)
          glob_fhist<<- rbind(glob_fhist,of)}
  }
  return(of)
}

# initialize global var for recording
init_glob_record_f <-function(){
  store_hist <<- TRUE # start history recording 
  nbcalls <<- 0
  glob_xhist<<- matrix(,1,ncol=nbvar)
  glob_fhist<<- matrix(,1,1)
}
