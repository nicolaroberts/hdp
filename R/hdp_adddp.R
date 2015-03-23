#' Add DPs to the HDP list object. 
#' 
#' @param hdp A HDP list object.
#' @param numdp The number of DPs to add (a positive integer). 
#' @param pp Integer specifying the index of the parental process for the new DPs.
#' @param cp Integer specifying the index of the concentration parameters for the new DPs.
#' @return A list with the updated HDP structure. 
#' @export
hdp_adddp <- function(hdp,numdp,pp,cp){
  
  HELDOUT <- 0L
  
  dpindex   <- hdp$numdp + 1:numdp
  hdp$numdp <- hdp$numdp + numdp
  hdp$dp <- c(hdp$dp, vector('list',numdp))
  for (ii in 1:numdp){
    jj <- dpindex[ii]
    tt <- hdp$conparam[[cp]]$numdp + 1
    hdp$dpstate[jj]              <- HELDOUT
    hdp$ppindex[jj]              <- as.integer(pp)
    hdp$cpindex[jj]              <- as.integer(cp)
    hdp$ttindex[jj]              <- as.integer(tt)
    hdp$conparam[[cp]]$numdp       <- as.integer(tt)
    hdp$conparam[[cp]]$totalnd[tt] <- 0L
    hdp$conparam[[cp]]$totalnt[tt] <- 0L
    hdp$dp[[jj]]$datacc            <- vector('integer')
    hdp$dp[[jj]]$classnd           <- 0L
    hdp$dp[[jj]]$classnt           <- 0L
    hdp$dp[[jj]]$beta              <- 1
    hdp$dp[[jj]]$alpha             <- vector('numeric')
    hdp$dp[[jj]]$numdata           <- 0L
    hdp$dp[[jj]]$datass            <- vector('integer')
  }
  
  return(hdp)
}
