#' Add concentration parameters to a hdpState object 
#' 
#' Add concentration parameters to a hdpState object by specifying the shape and
#' rate parameters of the gamma prior/s. DPs using these new concentration parameters
#' can be added with \code{\link{hdp_adddp}}. Data is assigned via \code{\link{hdp_setdata}}.
#' When initialised, the DP nodes are 'heldout' (not available for posterior sampling) 
#' and will need to be activated (see \code{\link{dp_activate}}). Finally, the posterior
#' sampling process (a Gibbs sampler) is run via \code{\link{hdp_posterior}}. 
#' 
#' @param hdp A hdpState object
#' @param alphaa Shape hyperparameters for the gamma priors over the DP concentration parameters.
#' @param alphab Rate hyperparameters for the gamma priors over the DP concentration parameters. 
#' @return A hdpState object updated with the new concentration parameters. See \code{\link{hdpState-class}} 
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_adddp}}, \code{\link{hdp_setdata}}, 
#'  \code{\link{dp_activate}}, \code{\link{hdp_posterior}}
#' @export
#' @examples
#' hdp_example <- hdp_init(c(0, 1, 1), c(1, 2, 2), rep(1, 6), rep(2, 2), rep(0.5, 2))
#' hdp_example <- hdp_addconparam(hdp_example, rep(1, 2), rep(1, 2))
#' 
hdp_addconparam <- function(hdp, alphaa, alphab){
  #input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  if (any(alphaa <= 0) | any(alphab <= 0)) {
    stop("alphaa and alphab must be positive")
  }
  if (length(alphaa) != length(alphab)) {
    stop("alphaa and alphab must have the same length")
  }
  
  # add new concentration parameters in
  old_numcp <- hdp@numconparam
  new_numcp <- length(alphaa)
  
  hdp@numconparam <- old_numcp + new_numcp
  hdp@conparam <- c(hdp@conparam, vector("list", new_numcp))
  
  # fill in conparam list
  for (cp in 1:new_numcp) {
    a <- alphaa[cp]
    b <- alphab[cp]
    numdpcp <- 0L
    hdp@conparam[[old_numcp+cp]] <- new("hdpConparam",
                              alphaa  = a,
                              alphab  = b,
                              numdp   = numdpcp,
                              alpha   = a / b,
                              totalnd = as.integer(rep(0, numdpcp)),
                              totalnt = as.integer(rep(0, numdpcp)))
  }
  
  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
