#' Freeze DP nodes
#' 
#' Freezes previously active DP nodes. A frozen DP node is not included in posterior sampling,
#' but its statistics \emph{are} considered in the sampling of other active DPs.
#' This is useful for conditioning on a previous dataset. First, set up a HDP
#' for one dataset, run the posterior sampling chain, and then freeze all old nodes
#' (except the top DP). Add new DP nodes with new data and run a
#' second posterior sampling chain over the new nodes (\emph{given} the information in the frozen nodes). 
#' 
#' @param hdp A hdpState object
#' @param dpindex Indices of the DPs to freeze
#' @return A hdpState object with the specified DP nodes frozen. See \code{\link{hdpState-class}} 
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_addconparam}}, \code{\link{hdp_adddp}}, 
#'  \code{\link{hdp_setdata}}, \code{\link{dp_activate}}, \code{\link{hdp_posterior}}
#' @export
dp_freeze <- function(hdp, dpindex){
  
  # input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  if (any(dpindex < 1) |
        any(dpindex > hdp@numdp) |
        any(dpindex %% 1 != 0) |
        any(duplicated(dpindex))) {
    stop("dpindex must be positive integers no greater than 
         numdp(hdp) with no duplicates")
  }
  
  ACTIVE  <- 2L
  FROZEN  <- 1L
  
  dpindex <- sort(dpindex)
  
  for (kk in 1:length(dpindex)){
    jj <- dpindex[kk] 
    if (hdp$dpstate[jj] != ACTIVE){
      stop("Can only freeze a DP that is activated")
    }

    hdp$dp[[jj]]$beta <- numeric(0)
    hdp$dp[[jj]]$alpha <- numeric(0)
    hdp$dpstate[jj]  <- FROZEN
  }
  
  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
