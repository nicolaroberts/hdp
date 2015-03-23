hdp_getstate <- function(hdp){
  hdpstate <- list()
  hdpstate$numclass <- hdp$base$numclass
  hdpstate$classqq  <- hdp$base$classqq
  hdpstate$classnd  <- vector('list',hdp$numdp)
  hdpstate$beta     <- vector('list',hdp$numdp)
  hdpstate$alpha    <- rep(0,hdp$numconparam)

  for (jj in 1:hdp$numdp){
    hdpstate$classnd[[jj]] <- hdp$dp[[jj]]$classnd
    hdpstate$beta[[jj]]    <- hdp$dp[[jj]]$beta
  }
  for (cp in 1:hdp$numconparam){
    hdpstate$alpha[cp]   <- hdp$conparam[[cp]]$alpha
  }
  return(hdpstate)
}


