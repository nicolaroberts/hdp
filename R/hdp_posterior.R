#' Posterior sampling chain across activated DPs.
#'
#' Run a Gibbs sampler over the activated DP nodes of a Hierarchichal Dirichlet Process.
#' Each iteration re-assigns the cluster allocation of every data item.
#' Run \code{burnin} iterations, and then collect \code{n} samples from the chain
#' with \code{space} iterations between each collected sample. To plot output,
#' see \code{\link{plot_lik}}, \code{\link{plot_numcluster}}, and
#' \code{\link{plot_data_assigned}}. Can collect multiple
#' independent HDP sampling chains in a hdpSampleMulti object via \code{\link{hdp_multi_chain}}.
#' Components are extracted via \code{\link{hdp_extract_components}}.
#'
#' @param hdp A hdpState object
#' @param burnin The number of burn-in iterations.
#' @param n The number of posterior samples to collect.
#' @param space The number of iterations between collected samples.
#' @param cpiter The number of iterations of concentration parameter sampling to perform after each iteration.
#' @param seed The (integer) seed that can be set to reproduce output. Default is a
#'  random seed from 1 -- 10^7, reported in the output.
#' @param verbosity Verbosity of debugging statements.
#'  0 (least verbose) -- 4 (most verbose). 0 highly recommended - only change for debugging small examples.
#' @return A hdpSampleChain object with the salient information from each
#'  posterior sample. See \code{\link{hdpSampleChain-class}}
#' @seealso \code{\link{hdp_multi_chain}}, \code{\link{hdp_extract_components}},
#'  \code{\link{cull_posterior_samples}}, \code{\link{plot_lik}}, \code{\link{plot_numcluster}},
#'  \code{\link{plot_data_assigned}}
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' my_hdp <- hdp_init(ppindex=0, cpindex=1, hh=rep(1, 6), alphaa=rep(1, 3), alphab=rep(2, 3))
#' my_hdp <- hdp_adddp(my_hdp, 2, 1, 2)
#' my_hdp <- hdp_adddp(my_hdp, 10, c(rep(2, 5), rep(3, 5)), 3)
#' my_hdp <- hdp_setdata(my_hdp, 4:13, example_data_hdp)
#' my_hdp <- dp_activate(my_hdp, 1:13, 2)
#' my_hdp_chain <- hdp_posterior(my_hdp, 100, 100, 10)
hdp_posterior <- function(hdp, burnin, n, space, cpiter=1,
                          seed=sample(1:10^7, 1), verbosity=0){

  # input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  for (arg in c("burnin", "n", "space", "cpiter")) {
    x <- get(arg)
    if (x < 1 | x %% 1 != 0) stop(paste(arg, "must be a positive integer"))
  }
  if (verbosity < 0 |
        verbosity > 4 |
        verbosity %% 1 != 0) stop("verbosity must be integer from 0--4")
  if (seed %% 1 != 0) stop("seed must be an integer")

  # set seed
  set.seed(seed, kind="Mersenne-Twister", normal.kind="Inversion")



  # function to return difference in time in minute units
  mindifftime <- function(t1, t2){
    as.numeric(t2-t1, units="mins")
  }

  # initialise likelihood vector
  totiter <- burnin + n * space
  lik     <- rep(0, totiter)



  # translate hdp hdpState (S4 class) to plain list so C code can parse
  hdplist <- as.list(hdp)

  # record start time
  starttime <- Sys.time()



  # run burn in iterations, update hdplist, fill in lik
  output <- iterate(hdplist, burnin, cpiter, verbosity)
  hdplist <- output[[1]]
  lik[1:burnin] <- output[[2]]



  #report burn-in time
  prevtime <- Sys.time()
  print(sprintf("%d burn-in iterations in %1.1f mins",
                burnin, mindifftime(starttime, prevtime)))
  curriter <- burnin

  # initialise list for posterior sample output
  sample  <- rep(list(hdp_getstate(hdplist)), n)

  # collect n posterior samples
  for (samp in 1:n){

    output <- iterate(hdplist, space, cpiter, verbosity)
    hdplist <- output[[1]]
    lik[burnin + (samp-1) * space + (1:space)] <- output[[2]]
    sample[[samp]] <- hdp_getstate(hdplist)

    #report time every 10 samples if > 1 min has passed
    tracktime <- Sys.time()
    curriter <- curriter + space
    if (mindifftime(prevtime, tracktime) > 1 & samp %% 10 == 0){
      elapsedtime <- mindifftime(starttime, tracktime)
      print(sprintf("time %1.1f ETC %1.1f mins",
                    elapsedtime, elapsedtime / curriter * totiter))
      prevtime <- tracktime
    }
  }



  numclass <- sapply(sample, function(x) x$numclass)
  classqq <- lapply(sample, function(x) x$classqq)
  classnd <- lapply(sample, function(x) as(x$classnd, "dgCMatrix"))
  alpha <- t(sapply(sample, function(x) x$alpha))

  # if only one conparam, then alpha can have wrong dims (vector not matrix)
  if (dim(alpha)[1]==1 & n > 1) {
    alpha <- matrix(alpha, ncol=1)
  }

  #translate hdplist back to HDPObject class
  hdp <- as.hdpState(hdplist)
  remove(hdplist)

  ans <- new("hdpSampleChain",
             seed = as.integer(seed),
             settings = list(burnin=burnin,
                             n=n,
                             space=space,
                             cpiter=cpiter),
             hdp = hdp,
             lik = lik,
             numcluster = numclass,
             cp_values = alpha,
             clust_categ_counts = classqq,
             clust_dp_counts = classnd,
             numcomp = as.integer(NULL),
             prop.ex = as.numeric(NULL),
             comp_cos_merge = as.numeric(NULL),
             comp_categ_counts = list(),
             comp_dp_counts = list(),
             comp_categ_distn = list(),
             comp_dp_distn = list())

  # check validity and return
  if (!validObject(ans)) warning("Not a valid hdpSampleChain object.")
  return(ans)
}
