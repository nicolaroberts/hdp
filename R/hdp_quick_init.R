#' Initialise a simple, default HDP structure
#'
#' Initialise a hdpState object with a basic default structure of one top parent DP node
#' with no associated data, and one child DP node per row of data. Every DP node
#' shares the same concentration parameter, and will automatically be 'activated'
#' (made available for posterior samplig). The base distribution is a uniform Dirichlet
#' with psuedocount 1 in each data category. Can immediately run \code{\link{hdp_posterior}}
#' to collect posterior samples.
#' To define a custom HDP structure, see \code{\link{hdp_init}} and \code{\link{hdp_prior_init}}.
#'
#' @param data A \code{data.frame} or \code{matrix} of counts with one row for every sample
#'  and one column for every data category.
#' @param initcc Number of initial data clusters
#'  (every data item is randomly assigned to a cluster to start with).
#' @param alphaa Shape hyperparameter for the gamma prior over the concentration parameter.
#' @param alphab Rate hyperparameter for the gamma prior over the concentration parameter.

#' @return A hdpState object with a basic default structure. See \code{\link{hdpState-class}}
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_posterior}}, \code{\link{hdp_prior_init}}
#' @export
#' @examples
#' my_quick_hdp <- hdp_quick_init(example_data_hdp)
#' my_quick_hdp_chain <- hdp_posterior(my_quick_hdp, 100, 50, 10, 5)

hdp_quick_init <- function(data, initcc=2, alphaa=1, alphab=1){

  # input checks
  if (!class(data) %in% c("matrix", "data.frame")) {
    stop("data must be data.frame or matrix")
  }
  if (any(data %% 1 != 0) | any(data < 0)) {
    stop("data must contain non-negative integer values")
  }
  if (initcc < 1 | initcc %% 1 != 0) stop("initcc must be a positive integer")
  if (alphaa <= 0 | alphab <= 0) {
    stop("alphaa and alphab must be positive")
  }
  if (length(alphaa) != 1 | length(alphab) != 1) {
    stop("alphaa and alphab must have length 1")
  }

  # number of data categories, samples, DPs
  numcateg <- ncol(data)
  nsamp <- as.integer(nrow(data))
  numdp <- as.integer(nsamp + 1)

  # hh base distribution is simple uniform Dirichlet
  hh <- rep(1, numcateg)

  # ppindex is one 0 then all 1
  ppindex <- c(0, rep(1, nsamp))

  # cpindex is all 1
  cpindex <- rep(1, numdp)

  # initialise hdpState
  hdp <- hdp_init(ppindex, cpindex, hh, alphaa, alphab)

  # add data
  hdp <- hdp_setdata(hdp, 2:numdp, data)

  # activate all DPs and initialise clusters
  hdp <- dp_activate(hdp, 1:numdp, initcc)

  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
