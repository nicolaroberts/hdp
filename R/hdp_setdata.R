#' Assign data to DPs in a HDP list object. 
#' The DPs having data assigned must still be 'held out'. 
#' @param hdp A HDP list object.
#' @param dpindex The indices of the DPs to assign data to.
#' @param data A data.frame with: one row for every sample to assign; one column for every category; counts of that data category in that sample.
#' @return A list with the updated HDP structure. 
#' @export
hdp_setdata <- function(hdp,dpindex,data){
  
  datass <- apply(data, 1, function(x) rep(1:ncol(data), x))

  check_dpindex(hdp$numdp,dpindex)
  if (length(dpindex) != length(datass)){
    stop('dpindex and datass lengths do not match')
  }
  
  HELDOUT <- 0L
  
  for (jj in 1:length(dpindex)){
    if (hdp$dpstate[dpindex[jj]] != HELDOUT){
      stop('Cannot set data for DPs that are not held out')
    }
    hdp$dp[[dpindex[jj]]]$numdata <- length(datass[[jj]])
    hdp$dp[[dpindex[jj]]]$datass  <- datass[[jj]]
  }
  
  return(hdp)
}

