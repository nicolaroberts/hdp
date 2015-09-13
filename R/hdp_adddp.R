#' Add DPs to a hdpState object 
#' 
#' Add DP nodes to a hdpState object and specify each parent relationship
#' and concentration parameter. Concentration parameters can be added to a 
#' hdpState object with \code{\link{hdp_addconparam}}. Data is assigned via \code{\link{hdp_setdata}}.
#' When initialised, the DP nodes are 'heldout' (not available for posterior sampling) 
#' and will need to be activated (see \code{\link{dp_activate}}). Finally, the posterior
#' sampling process (a Gibbs sampler) is run via \code{\link{hdp_posterior}}. 
#' 
#' 
#' @param hdp A hdpState object 
#' @param numdp The number of DPs to add 
#' @param ppindex Index (or indices) of the parental process(es) for the new DPs.
#' @param cpindex Index (or indices) of the concentration parameters for the new DPs.
#' @return A hdpState object with the updated HDP structure. See \code{\link{hdpState-class}}
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_setdata}}, 
#'  \code{\link{dp_activate}}, \code{\link{hdp_posterior}}
#' @export
#' @examples
#' my_hdp <- hdp_init(ppindex=0, cpindex=1, hh=rep(1, 6), alphaa=rep(1, 3), alphab=rep(2, 3))
#' # add two more DPs with parent '1' and concentration parameter '2'
#' my_hdp <- hdp_adddp(my_hdp, 2, 1, 2)
#' my_hdp
#' 
#' hdp_example <- hdp_init(c(0, 1, 1), c(1, 2, 2), rep(1, 6), rep(2, 2), rep(0.5, 2))
#' # add six more DPs, three with parent '2', three with parent '3', 
#' # and all with concentration parameter '2'
#' hdp_example <- hdp_adddp(hdp_example, 6, c(2, 2, 2, 3, 3, 3), 2)
#' hdp_example
hdp_adddp <- function(hdp, numdp, ppindex, cpindex){
  
  # input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  if (numdp < 1 | numdp %% 1 != 0) stop("numdp must be positive integer")
  # adjust length of ppindex and cpindex if single integer
  if (length(ppindex)==1 & numdp>1) ppindex <- rep(ppindex, numdp)
  if (length(cpindex)==1 & numdp>1) cpindex <- rep(cpindex, numdp)
  if (length(ppindex) != numdp | length(cpindex) != numdp){
    stop("ppindex and cpindex must be either a single integer,
         or have length equal to numdp")
  }
  if (any(ppindex < 0) |
        any(ppindex %% 1 != 0) |
        any(ppindex >= (hdp@numdp+1):(hdp@numdp+numdp))) {
    stop("ppindex must be non-negative integer/s,
         referring to a parent of smaller index")
  }
  if (any(cpindex < 1) | 
        any(cpindex %% 1 != 0) | 
        any(cpindex > hdp@numconparam) |
        length(cpindex) != length(ppindex)) {
    stop("cpindex must be positive integer/s, no greater than the number of 
         concentration parameters, and same length as ppindex")
  } 
  
  
  HELDOUT <- 0L
  
  # add new DPs
  dpindex   <- hdp@numdp + 1:numdp
  hdp@numdp <- hdp@numdp + as.integer(numdp)
  hdp@dp <- c(hdp@dp, vector("list", numdp))
  for (ii in 1:numdp){
    jj <- dpindex[ii]
    pp <- ppindex[ii]
    cp <- cpindex[ii]
    tt <- hdp@conparam[[cp]]@numdp + 1
    hdp@dpstate[jj]              <- HELDOUT
    hdp@ppindex[jj]              <- as.integer(pp)
    hdp@cpindex[jj]              <- as.integer(cp)
    hdp@ttindex[jj]              <- as.integer(tt)
    hdp@conparam[[cp]]@numdp       <- as.integer(tt)
    hdp@conparam[[cp]]@totalnd[tt] <- 0L
    hdp@conparam[[cp]]@totalnt[tt] <- 0L
    hdp@dp[[jj]] <- new("hdpDP",
                        datacc  = vector("integer"),
                        classnd = 0L,
                        classnt = 0L,
                        beta    = 1,
                        alpha   = vector("numeric"),
                        numdata = 0L,
                        datass  = vector("integer"))
  }
  
  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
