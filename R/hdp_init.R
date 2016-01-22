#' Initialise a hdpState object
#'
#' Initialise a hdpState object with one or more DP nodes and their parent relationships,
#' the parameters of the base Dirichlet distribution, and a set of hyperparameters
#' for the gamma priors over the DP concentration parameters. Further DP nodes can
#' be added with \code{\link{hdp_adddp}}, and further concentration parameters can
#' be added with \code{\link{hdp_addconparam}}. Data is assigned via \code{\link{hdp_setdata}}.
#' When initialised, the DP nodes are 'heldout' (not available for posterior sampling)
#' and will need to be activated (see \code{\link{dp_activate}}). Finally, the posterior
#' sampling process (a Gibbs sampler) is run via \code{\link{hdp_posterior}}.
#'
#'
#' @param ppindex Index (or indices) of the parental process(es) for the initial DPs.
#'  The 'top' DP should have parent process '0' (the base Dirichlet distribution).
#' @param cpindex Index (or indices) of the concentration parameter(s) for the initial DPs.
#' @param hh Parameters of the base Dirichlet distribution (like psuedocounts across categories).
#'  Must be a vector with length equal to the number of data item categories.
#' @param alphaa Shape hyperparameters for the gamma priors over the DP concentration parameters.
#' @param alphab Rate hyperparameters for the gamma priors over the DP concentration parameters.
#' @return A hdpState object with the initial HDP structure. See \code{\link{hdpState-class}}
#' @seealso \code{\link{hdp_quick_init}}, \code{\link{hdp_prior_init}},
#'  \code{\link{hdp_addconparam}}, \code{\link{hdp_adddp}},
#'  \code{\link{hdp_setdata}}, \code{\link{dp_activate}}, \code{\link{hdp_posterior}}
#' @export
#' @examples
#' # initialise a HDP with just one 'top' DP node off the base distribution,
#' # a uniform Dirichlet base distribution over six possible data categories,
#' # and three possible concentration parameters to be shared across the HDP tree
#' # (top DP using conparam number 1), each with hyperparameters (1,2).
#' hdp_init(ppindex=0, cpindex=1, hh=rep(1, 6), alphaa=rep(1, 3), alphab=rep(2, 3))
#'
#' # initialise a HDP with one 'top' DP node off the base distribution,
#' # AND two children DP nodes off that parent. The two children DPs share a different
#' # concentration parameter (hyperparameters are (2, 0.5)).
#' hdp_init(ppindex=c(0, 1, 1), cpindex=c(1, 2, 2), hh=rep(1, 6), alphaa=rep(2, 2), alphab=rep(0.5, 2))

hdp_init <- function(ppindex, cpindex, hh, alphaa, alphab){

  # input checks
  if (any(ppindex < 0) |
        any(ppindex %% 1 != 0) |
        any(ppindex >= 1:length(ppindex))) {
    stop("ppindex must be non-negative integer/s,
         referring to a parent of smaller index")
  }
  if (any(cpindex < 1) |
        any(cpindex %% 1 != 0) |
        any(cpindex > length(alphaa)) |
        length(cpindex) != length(ppindex)) {
    stop("cpindex must be positive integer/s, no greater than
         the length of alphaa and alphab, and same length as ppindex")
  }
  if (any(hh <= 0)) stop("hh must be a positive vector")
  if (any(alphaa <= 0) | any(alphab <= 0)) {
    stop("alphaa and alphab must be positive")
  }
  if (length(alphaa) != length(alphab)) {
    stop("alphaa and alphab must have the same length")
  }

  HELDOUT <- 0

  # initialise base distribution
  base <- new("hdpBase",
              hh = hh,
              classqq = newclass(hh),
              numclass = 0L)

  # initialise hdpState
  numdp <- length(ppindex)
  numconparam <- length(alphaa)
  hdp <- new("hdpState",
             numdp         = numdp,
             numconparam   = numconparam,
             base          = base,
             conparam      = vector("list", numconparam),
             dp            = vector("list", numdp),
             dpstate       = as.integer(HELDOUT * rep(1, numdp)),
             ppindex       = as.integer(ppindex),
             cpindex       = as.integer(cpindex),
             ttindex       = as.integer(rep(0,numdp)),
             initcc        = as.integer(NULL),
             seed_activate = as.integer(NULL),
             pseudoDP      = as.integer(NULL)
             )

  # fill in ttindex slot
  tt <- rep(0, hdp@numconparam)
  for (jj in 1:hdp@numdp){
    cp              <- cpindex[jj]
    tt[cp]          <- tt[cp] + 1
    hdp@ttindex[jj] <- as.integer(tt[cp])
  }

  # fill in dp list
  for (jj in 1:hdp@numdp){
    hdp@dp[[jj]] <- new("hdpDP",
                        datacc  = vector("integer"),
                        classnd = 0L,
                        classnt = 0L,
                        beta    = 1,
                        alpha   = vector("numeric"),
                        numdata = 0L,
                        datass  = vector("integer"))
  }

  # fill in conparam list
  for (cp in 1:hdp@numconparam){
    a <- alphaa[cp]
    b <- alphab[cp]
    numdpcp <- as.integer(sum(cpindex == cp))
    hdp@conparam[[cp]] <- new("hdpConparam",
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
