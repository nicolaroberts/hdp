#' Initialise a HDP structure incorporating prior knowledge
#'
#' Initialise a hdpState object incorporating prior knowlegde of some components
#'(categorical data distributions). The structure has one top parent DP node
#' with no associated data ('active' and available for posterior sampling),
#' and one child DP node per prior component ('frozen' and held out
#' from posterior sampling).
#'
#' @param prior_distn Matrix of prior distributions (columns must each sum to 1,
#'  number of rows matches number of data categories)
#' @param prior_pseudoc Vector of pseudocounts contributed by each prior distribution
#' @param hh Parameters of the base Dirichlet distribution. Must be a vector with length equal
#'  to the number of data item categories.
#' @param alphaa Shape hyperparameters for the gamma priors over the DP concentration parameters.
#' @param alphab Rate hyperparameters for the gamma priors over the DP concentration parameters.
#' @seealso \code{\link{hdp_init}}
#' @return A hdpState object with one frozen node per prior component. See \code{\link{hdpState-class}}
#' @export
#' @examples
#' # example dataset with 10 data categories, and 100 samples.
#' # Two components are known a priori.
#' hdp <- hdp_prior_init(example_known_priors, rep(1000, 2), hh=rep(1, 10),
#'              alphaa=c(1,1), alphab=c(1,1))
#' hdp <- hdp_addconparam(hdp, alphaa=c(1,1), alphab=c(1,1))
#' hdp <- hdp_adddp(hdp, 101, c(1, rep(4, 100)), c(3, rep(4, 100)))
#' hdp <- hdp_setdata(hdp, 5:104, example_data_hdp_prior)
#' hdp <- dp_activate(hdp, 4:104, initcc=4, seed=81479)
#' hdp <- hdp_posterior(hdp, burnin=2000, n=50, space=50, cpiter=3, seed=1e6)
#' hdp_ex <- hdp_extract_components(hdp)
#' plot_comp_size(hdp_ex)
#' plot_comp_distn(hdp_ex)
#' plot_dp_comp_exposure(hdp_ex, 5:104, col_comp=rainbow(5))


hdp_prior_init <- function(prior_distn, prior_pseudoc, hh, alphaa, alphab){
  # check input
  if (!is.matrix(prior_distn) | any(prior_distn<0) |
      any(abs(colSums(prior_distn)-1) > 1e-3)) {
    stop("prior_distn must be a matrix of probability vectors, with columns summing to 1")
  }
  if (!is.vector(prior_pseudoc) | any(prior_pseudoc<0) |
      length(prior_pseudoc)!=ncol(prior_distn)) {
    stop("prior_pseudoc must be a non-negative vector with length equal to ncol(prior_distn)")
  }
  if (!is.vector(hh) | any(hh <= 0) | length(hh)!=nrow(prior_distn)) {
    stop("hh must be a positive vector with length equal to nrow(prior_distn)")
  }
  if (any(alphaa <= 0) | any(alphab <= 0)) {
    stop("alphaa and alphab must be positive")
  }
  if (length(alphaa) != 2 | length(alphab) != 2) {
    stop("alphaa and alphab must have length 2")
  }

  FROZEN <- 1L
  ACTIVE <- 2L

  # initialse HDP structure and assign fake pseudodata
  nsig <- ncol(prior_distn)
  hdp <- hdp_init(ppindex=c(0, rep(1, nsig)),
                  cpindex=c(1, rep(2, nsig)),
                  hh=hh,
                  alphaa=alphaa,
                  alphab=alphab)
  hdp@pseudoDP <- 1L + (1:nsig)

  prior_mut_count <- t(round(prior_distn %*% diag(prior_pseudoc)))
  hdp <- hdp_setdata(hdp, 2:numdp(hdp), prior_mut_count)

  # initialize numclass and classqq
  hdp <- qq_addclass(hdp, nsig)

  # initialize clustering of HDP with each node of fake pseudodata frozen as a fixed cluster
  for (jj in 1:numdp(hdp)){
    pp <- hdp@ppindex[jj]
    cp <- hdp@cpindex[jj]
    tt <- hdp@ttindex[jj]

    hdp@dp[[jj]]@datacc <- as.integer(rep(jj-1, hdp@dp[[jj]]@numdata))
    hdp@base@classqq <- additems(hdp@base@classqq, hdp@dp[[jj]]@datacc,
                                 hdp@dp[[jj]]@datass)
    hdp@dp[[jj]]@classnd <- tabulate(hdp@dp[[jj]]@datacc,
                                     nbins=hdp@base@numclass+1)
    hdp@dp[[jj]]@classnt <- 0L


    if (pp == 0){
      hdp@dp[[jj]]@beta <- rep(1,hdp@base@numclass+1) / (hdp@base@numclass+1)
    } else {
      hdp@dp[[jj]]@beta <- hdp@dp[[pp]]@beta
    }
    hdp@dp[[jj]]@alpha <- hdp@conparam[[cp]]@alpha
  }

  # calculate tables, and FREEZE pseudodata nodes
  for (jj in numdp(hdp):1){
    pp <- hdp@ppindex[jj]
    cp <- hdp@cpindex[jj]
    tt <- hdp@ttindex[jj]
    alpha <- hdp@dp[[jj]]@alpha

    if (pp == 0){
      hdp@dp[[jj]]@classnt <- as.integer(hdp@dp[[jj]]@classnd > 0)
      hdp@dpstate[jj]  <- ACTIVE
    } else {
      hdp@dp[[pp]]@classnd <- as.integer(hdp@dp[[pp]]@classnd -
                                           hdp@dp[[jj]]@classnt)
      hdp@dp[[jj]]@classnt <- as.integer(randnumtable(alpha*hdp@dp[[pp]]@beta,
                                                      hdp@dp[[jj]]@classnd))
      hdp@dp[[pp]]@classnd <- as.integer(hdp@dp[[pp]]@classnd +
                                           hdp@dp[[jj]]@classnt)

      hdp@dpstate[jj]  <- FROZEN
    }

    hdp@conparam[[cp]]@totalnd[tt] <- as.integer(sum(hdp@dp[[jj]]@classnd))
    hdp@conparam[[cp]]@totalnt[tt] <- as.integer(sum(hdp@dp[[jj]]@classnt))
  }

  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}

