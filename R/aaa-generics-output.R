# Generics for hdpSampleChain class -------------

#' @describeIn hdpSampleChain Convert to list class
#' @export
#' @param x Object of class hdpSampleChain
#' @param ... unused
setMethod("as.list",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- list(seed=x@seed,
                        settings=x@settings,
                        hdp=x@hdp,
                        lik=x@lik,
                        numcluster=x@numcluster,
                        cp_values=x@cp_values,
                        clust_categ_counts=x@clust_categ_counts,
                        clust_dp_counts=x@clust_dp_counts,
                        numcomp=x@numcomp,
                        prop.ex=x@prop.ex,
                        comp_cos_merge=x@comp_cos_merge,
                        comp_categ_counts=x@comp_categ_counts,
                        comp_dp_counts=x@comp_dp_counts,
                        comp_categ_distn=x@comp_categ_distn,
                        comp_dp_distn=x@comp_dp_distn)
            return(ans)
          })


# show method
setMethod("show",
          "hdpSampleChain",
          function(object) {
            cat("Object of class", class(object), "\n")
            cat(" Random seed:", object@seed, "\n")
            cat(" Burn-in iters:", object@settings$burnin, "\n")
            cat(" Posterior samples:", object@settings$n, "\n")
            cat(" Iters b/w samples:", object@settings$space, "\n")
            cat(" Tabulate raw cluster number across samples:")
            print(table(object@numcluster))
            if (length(object@numcomp) != 1) {
              cat(" Components: NO. Run hdp_extract_components \n")
            } else {
              cat(" Components: YES. Prop of data explained =",
                  object@prop.ex,
                  " Merged if cosine sim >",
                  object@comp_cos_merge, "\n")
              cat(" Component number:", object@numcomp, "\n")
            }
            cat(" ----------\n")
            cat(" Final hdpState: \n")
            print(object@hdp)
          })

# Generics for hdpSampleMulti class -------------

#' @describeIn hdpSampleMulti Convert to list class
#' @export
#' @param x Object of class hdpSampleMulti
#' @param ... unused
setMethod("as.list",
          signature = "hdpSampleMulti",
          definition = function(x, ...) {
            ans <- list(chains=x@chains,
                        numcomp=x@numcomp,
                        prop.ex=x@prop.ex,
                        comp_cos_merge=x@comp_cos_merge,
                        comp_categ_counts=x@comp_categ_counts,
                        comp_dp_counts=x@comp_dp_counts,
                        comp_categ_distn=x@comp_categ_distn,
                        comp_dp_distn=x@comp_dp_distn)
            return(ans)
          })


# show method
setMethod("show",
          "hdpSampleMulti",
          function(object) {

            totalsamp <- sum(sapply(object@chains, function(x) x@settings$n))

            cat("Object of class", class(object), "\n")
            cat(" Number of chains:", length(object@chains), "\n")
            cat(" Total posterior samples:", totalsamp, "\n")
            if (length(object@numcomp) != 1) {
              cat(" Components: NO. Run hdp_extract_components \n")
            } else {
              cat(" Components: YES. Prop of data explained =",
                  object@prop.ex,
                  " Merged if cosine sim >",
                  object@comp_cos_merge, "\n")
              cat(" Component number:", object@numcomp, "\n")
            }
            cat(" ----------\n")
            cat(" Final hdpState from first chain: \n")
            print(object@chains[[1]]@hdp)
          })

# getter methods for hdpSampleChain and hdpSampleMulti objects -----------

setGeneric("sampling_seed", function(x, ...) standardGeneric("sampling_seed"))
#' @describeIn hdpSampleChain Get random seed used by \code{hdp_posterior}
#' @aliases sampling_seed
#' @export
setMethod("sampling_seed",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@seed
            return(ans)
          })

setGeneric("hdp_settings", function(x, ...) standardGeneric("hdp_settings"))
#' @describeIn hdpSampleChain Get settings of posterior sampling chain
#' @aliases hdp_settings
#' @export
setMethod("hdp_settings",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@settings
            return(ans)
          })

setGeneric("final_hdpState", function(x, ...) standardGeneric("final_hdpState"))
#' @describeIn hdpSampleChain Get hdpState object from the end of the posterior sampling chain
#' @aliases final_hdpState
#' @export
setMethod("final_hdpState",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@hdp
            return(ans)
          })

setGeneric("chains", function(x, ...) standardGeneric("chains"))
#' @describeIn hdpSampleMulti Get list of hdpSampleChain objects
#' @aliases chains
#' @export
setMethod("chains",
          signature = "hdpSampleMulti",
          definition = function(x, ...) {
            ans <- x@chains
            return(ans)
          })

setGeneric("lik", function(x, ...) standardGeneric("lik"))
#' @describeIn hdpSampleChain Get likelihood of data given model over all iterations
#' @aliases lik
#' @export
setMethod("lik",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@lik
            return(ans)
          })

setGeneric("numcluster", function(x, ...) standardGeneric("numcluster"))
#' @describeIn hdpSampleChain Get the number of clusters for each posterior sample
#' @aliases numcluster
#' @export
setMethod("numcluster",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@numcluster
            return(ans)
          })

setGeneric("cp_values", function(x, ...) standardGeneric("cp_values"))
#' @describeIn hdpSampleChain Get matrix of concentration parameter values for each posterior sample
#' @aliases cp_values
#' @export
setMethod("cp_values",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@cp_values
            return(ans)
          })

setGeneric("clust_categ_counts",
           function(x, ...) standardGeneric("clust_categ_counts"))
#' @describeIn hdpSampleChain Get category vs cluster counts for each posterior sample
#' @aliases clust_categ_counts
#' @export
setMethod("clust_categ_counts",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@clust_categ_counts
            return(ans)
          })

setGeneric("clust_dp_counts",
           function(x, ...) standardGeneric("clust_dp_counts"))
#' @describeIn hdpSampleChain Get dp node vs cluster counts for each posterior sample
#' @aliases clust_dp_counts
#' @export
setMethod("clust_dp_counts",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@clust_dp_counts
            return(ans)
          })


#' Get number of extracted components
#' @param x hdpSampleChain or hdpSampleMulti
#' @return number of components
#' @aliases numcomp
#' @export
setGeneric("numcomp",
           function(x) standardGeneric("numcomp"))

#' @describeIn hdpSampleChain Get number of extracted components for hdpSampleChain
setMethod("numcomp",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@numcomp
            return(ans)
          })

#' @describeIn hdpSampleMulti Get number of extracted components for hdpSampleMulti
setMethod("numcomp",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@numcomp
            return(ans)
          })

#' Get proportion of dataset explained (on average)
#' @param x hdpSampleChain or hdpSampleMulti
#' @return number of components
#' @aliases prop.ex
#' @export
setGeneric("prop.ex",
           function(x) standardGeneric("prop.ex"))

#' @describeIn hdpSampleChain Get proportion of dataset explained (on average) for hdpSampleChain
setMethod("prop.ex",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@prop.ex
            return(ans)
          })

#' @describeIn hdpSampleMulti Get proportion of dataset explained (on average) for hdpSampleMulti
setMethod("prop.ex",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@prop.ex
            return(ans)
          })

#' Get cos.merge setting
#' @param x hdpSampleChain or hdpSampleMulti
#' @return number of components
#' @aliases comp_cos_merge
#' @export
setGeneric("comp_cos_merge",
           function(x) standardGeneric("comp_cos_merge"))

#' @describeIn hdpSampleChain Get cos.merge setting for hdpSampleChain
setMethod("comp_cos_merge",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@comp_cos_merge
            return(ans)
          })

#' @describeIn hdpSampleMulti Get cos.merge setting for hdpSampleMulti
setMethod("comp_cos_merge",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@comp_cos_merge
            return(ans)
          })


#' Get sample vs category counts for each component
#' @param x hdpSampleChain or hdpSampleMulti
#' @return List of matrices (one for each component)
#'  counting the sample-category data assignment across all DP nodes.
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of data categories.
#' @aliases comp_categ_counts
#' @export
setGeneric("comp_categ_counts",
           function(x) standardGeneric("comp_categ_counts"))

#' @describeIn hdpSampleChain Get sample vs category counts for each component
setMethod("comp_categ_counts",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@comp_categ_counts
            return(ans)
          })

#' @describeIn hdpSampleMulti Get sample vs category counts for each component
setMethod("comp_categ_counts",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@comp_categ_counts
            return(ans)
          })


#' Get sample vs component counts for each DP
#' @param x hdpSampleChain or hdpSampleMulti
#' @return List of matrices (one for each DP)
#'  counting sample-component assignment (aggregating across data categories).
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of components.
#' @aliases comp_dp_counts
#' @export
setGeneric("comp_dp_counts",
           function(x) standardGeneric("comp_dp_counts"))

#' @describeIn hdpSampleChain Get sample vs component counts for each DP
setMethod("comp_dp_counts",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@comp_dp_counts
            return(ans)
          })

#' @describeIn hdpSampleMulti Get sample vs component counts for each DP
setMethod("comp_dp_counts",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@comp_dp_counts
            return(ans)
          })

#' Get mean distribution over data categories for each component
#' @param x hdpSampleChain or hdpSampleMulti
#' @return List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95\% credibility interval) distribution
#'  over data categories for each component. Number of rows is the number of
#'  components, and number of columns is the number of data categories. Rows sum to 1.
#' @aliases comp_categ_distn
#' @export
setGeneric("comp_categ_distn",
           function(x) standardGeneric("comp_categ_distn"))

#' @describeIn hdpSampleChain Get mean distribution over data categories for each component
setMethod("comp_categ_distn",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@comp_categ_distn
            return(ans)
          })

#' @describeIn hdpSampleMulti Get mean distribution over data categories for each component
setMethod("comp_categ_distn",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@comp_categ_distn
            return(ans)
          })


#' Get mean distribution over components for each DP
#' @param x hdpSampleChain or hdpSampleMulti
#' @return List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95\% credibility interval) distribution
#'  over components for each DP. Number of rows is the number of
#'  DPs, and number of columns is the number of components. Rows sum to 1.
#' @aliases comp_dp_distn
#' @export
setGeneric("comp_dp_distn",
           function(x) standardGeneric("comp_dp_distn"))

#' @describeIn hdpSampleChain Get mean distribution over components for each DP
setMethod("comp_dp_distn",
          signature = "hdpSampleChain",
          definition = function(x) {
            ans <- x@comp_dp_distn
            return(ans)
          })

#' @describeIn hdpSampleMulti Get mean distribution over components for each DP
setMethod("comp_dp_distn",
          signature = "hdpSampleMulti",
          definition = function(x) {
            ans <- x@comp_dp_distn
            return(ans)
          })

# Component extraction --------------

#' Extract major components from the raw clusters
#'
#' If prior components included via \code{\link{hdp_prior_init}}, they will be
#' preserved by \code{hdp_extract_components} and prefixed with "P".
#' Any new components in this case are prefixed with "N".
#'
#' @param x hdpSampleChain or hdpSampleMulti object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90)
#' @param redo Logical. If true - constituent chains with previously calculated components
#'  will be re-calculated. Only used for hdpSampleMulti.
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases hdp_extract_components
#' @seealso \code{\link{hdp_posterior}}, \code{\link{hdp_multi_chain}},
#'  \code{\link{plot_comp_size}}, \code{\link{plot_comp_distn}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @import clue
#' @include hdp_extract_comp_single.R hdp_extract_comp_multi.R
#' @export
#' @examples
#' hdp_extract_components(mut_example_chain)

setGeneric("hdp_extract_components",
           function(x, cos.merge=0.9, redo=TRUE)
             standardGeneric("hdp_extract_components"))

#' @describeIn hdp_extract_components Extract components for hdpSampleChain
setMethod("hdp_extract_components",
          signature = "hdpSampleChain",
          definition = function(x, cos.merge) {
            ans <- hdp_extract_comp_single(x, cos.merge)
            return(ans)
          })

#' @describeIn hdp_extract_components Extract components for hdpSampleMulti
setMethod("hdp_extract_components",
          signature = "hdpSampleMulti",
          definition = function(x, cos.merge, redo) {
            ans <- hdp_extract_comp_multi(x, cos.merge, redo)
            return(ans)
          })
