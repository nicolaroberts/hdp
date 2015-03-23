# posterior MCMC sampling for hdp data struct
# returns an R list with two elements: 
# the updated hdp data struct, and
# the vector of likelihoods for these iterations

#' @useDynLib hdp hdpMultinomial_iterate
iterate <- function(hdp,numiter,doconparam,dolik,dodebug){
  out <- .Call(hdpMultinomial_iterate, hdp, numiter, doconparam, dolik, dodebug, PACKAGE='hdp')
  return(out)
}