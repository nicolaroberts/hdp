#' Cull early posterior samples from a hdpSampleChain object
#'
#' Extend the 'burn-in' period and reduce the number of posterior samples taken
#' from a sampling chain by culling the first \code{ncull} posterior samples.
#' If components have been previously calculated for this sampling chain, they
#' will be removed and must be recalculated.
#'
#'
#' @param chain A hdpSampleChain object
#' @param ncull The number of posterior samples to cull
#' @return A hdpSampleChain object with the designated 'burn-in' period extended,
#'  and the number of posterior samples reduced by \code{ncull}
#' @export
#' @seealso \code{\link{plot_lik}}, \code{\link{plot_numcluster}},
#'  \code{\link{plot_data_assigned}}
#' @examples
#' chain_adj <- cull_posterior_samples(mut_example_chain, 20)
#' plot_lik(chain_adj)
#' plot_numcluster(chain_adj)
#' plot_data_assigned(chain_adj)
cull_posterior_samples <- function(chain, ncull){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (ncull < 1 | ncull >= hdp_settings(chain)$n |
        ncull %% 1 != 0) {
    stop("ncull must be positive integer less than number of samples in chain")
  }

  # if component slot exist, clear them
  if (length(comp_categ_counts(chain)) > 0 ){
    message("Deleting component calculations. Re-run hdp_extract_components")
    chain@comp_settings <- list()
    chain@comp_categ_counts <- list()
    chain@comp_dp_counts <- list()
  }

  old_n <- chain@settings$n
  old_burnin <- chain@settings$burnin
  space <- chain@settings$space

  # adjust settings
  chain@settings$n <- old_n - ncull
  chain@settings$burnin <- old_burnin + ncull*space

  # cut first ncull elements of numcluster, cp_values,
  # clust_categ_counts, clust_dp_counts.

  chain@numcluster <- chain@numcluster[-c(1:ncull)]
  chain@cp_values <- as.matrix(chain@cp_values[-c(1:ncull),])
  chain@clust_categ_counts <- chain@clust_categ_counts[-c(1:ncull)]
  chain@clust_dp_counts <- chain@clust_dp_counts[-c(1:ncull)]

  # check validity and return
  if (!validObject(chain)) warning("Not a valid hdpSampleChain object.")
  return(chain)

}
