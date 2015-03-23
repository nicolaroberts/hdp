#' Freezes previously active DPs. A frozen DP is not included in posterior sampling, but the statistics for this DP are used in the sampling of other active DPs.
#' @param hdp A HDP list object. 
#' @param dpindex A vector of indices indicating which DPs in the HDP object to freeze. 
#' @export
dp_freeze <- function(hdp,dpindex){
  
  check_dpindex(hdp$numdp,dpindex)
  
  ACTIVE  <- 2L
  FROZEN  <- 1L
  
  dpindex <- sort(dpindex)
  # initialize state of HDP
  for (kk in 1:length(dpindex)){
    jj <- dpindex[kk] 
    if (hdp$dpstate[jj] != ACTIVE){
      stop('Can only freeze a DP that is activated')
    }

    hdp$dp[[jj]]$beta <- numeric(0)
    hdp$dp[[jj]]$alpha <- numeric(0)
    hdp$dpstate[jj]  <- FROZEN
  }
  return(hdp)
}


