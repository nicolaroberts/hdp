# extract key info from this hdp iteration (list class)
# numclass = number of clusters
# classqq = category vs cluster counts overall (matrix)
# classnd = dp vs cluster counts (matrix)
# beta = dp vs cluster weights (matrix)
# alpha = conparam values (vector)

hdp_getstate <- function(hdp){
  hdpstate <- list()
  hdpstate$numclass <- hdp$base$numclass
  hdpstate$classqq  <- hdp$base$classqq
  hdpstate$classnd  <- t(sapply(hdp$dp, function(x) x$classnd))
  hdpstate$beta     <- t(sapply(hdp$dp, function(x) x$beta))
  hdpstate$alpha    <- sapply(hdp$conparam, function(x) x$alpha)
  return(hdpstate)
}
