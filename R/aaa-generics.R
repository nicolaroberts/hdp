# Generics for hdpBase class -------------

#' @describeIn hdpBase convert to list class
#' @param x Object of class hdpBase
#' @param ... unused
#' @export
setMethod("as.list",
          signature = "hdpBase",
          definition = function(x, ...) {
            ans <- list(hh=x@hh,
                        classqq=x@classqq,
                        numclass=x@numclass)
            return(ans)
          })


# convert list to hdpBase class
setGeneric("as.hdpBase", function(x, ...) standardGeneric("as.hdpBase"))
setMethod("as.hdpBase",
          signature = "list",
          definition = function(x, ...) {
            ans <- new("hdpBase",
                       hh=x$hh,
                       classqq=x$classqq,
                       numclass=x$numclass)
            return(ans)
          })


# show method
setMethod("show",
          "hdpBase",
          function(object) {
            numcat <- length(object@hh)
            numclass <- object@numclass
            cat("Object of class", class(object), "\n")
            cat(" Number of data categories:", numcat, "\n")
            if (numcat <= 5){
              cat(" Prior over data categories:", object@hh, "\n")
            } else {
              cat(" Prior over data categories:", object@hh[1:5], "...\n")
            }
            cat(" Number of clusters:", numclass, "\n")
          })


# Generics for hdpConparam class -------------

#' @describeIn hdpConparam convert to list class
#' @export
#' @param x Object of class hdpConparam
#' @param ... unused
setMethod("as.list",
          signature = "hdpConparam",
          definition = function(x, ...) {
            ans <- list(alphaa = x@alphaa,
                        alphab = x@alphab,
                        numdp = x@numdp,
                        alpha = x@alpha,
                        totalnd = x@totalnd,
                        totalnt = x@totalnt)
            return(ans)
          })


# convert list to hdpConparam class
setGeneric("as.hdpConparam", function(x, ...) standardGeneric("as.hdpConparam"))
setMethod("as.hdpConparam",
          signature = "list",
          definition = function(x, ...) {
            ans <- new("hdpConparam",
                       alphaa = x$alphaa,
                       alphab = x$alphab,
                       numdp = x$numdp,
                       alpha = x$alpha,
                       totalnd = x$totalnd,
                       totalnt = x$totalnt)
            return(ans)
          })


# show method
setMethod("show",
          "hdpConparam",
          function(object) {
            cat("Object of class", class(object), "\n")
            cat(" Parameters of gamma prior:",
                object@alphaa, object@alphab, "\n")
            cat(" Value:", object@alpha, "\n")
            cat(" Number of DPs using this concentration param:",
                object@numdp, "\n")
          })


# Generics for hdpDP class -------------

#' @describeIn hdpDP convert to list class
#' @export
#' @param x Object of class hdpDP
#' @param ... unused
setMethod("as.list",
          signature = "hdpDP",
          definition = function(x, ...) {
            ans <- list(datacc=x@datacc,
                        classnd=x@classnd,
                        classnt=x@classnt,
                        beta=x@beta,
                        alpha=x@alpha,
                        numdata=x@numdata,
                        datass=x@datass)
            return(ans)
          })


# convert list to hdpDP class
setGeneric("as.hdpDP", function(x, ...) standardGeneric("as.hdpDP"))
setMethod("as.hdpDP",
          signature = "list",
          definition = function(x, ...) {
            ans <- new("hdpDP",
                       datacc=x$datacc,
                       classnd=x$classnd,
                       classnt=x$classnt,
                       beta=x$beta,
                       alpha=x$alpha,
                       numdata=x$numdata,
                       datass=x$datass)
            return(ans)
          })


# show method
setMethod("show",
          "hdpDP",
          function(object) {
            cat("Object of class", class(object), "\n")
            cat(" Number of data items:", object@numdata, "\n")
            if (object@numdata==0){
              cat(" Number of child data tables in each cluster:",
                  object@classnd, "\n")
            } else {
              cat(" Number of data items in each cluster:",
                  object@classnd, "\n")
              cat(" Number of data items in each category:",
                  table(object@datass), "\n")
            }
          })


# Generics for hdpState class -------------

#' @describeIn hdpState Convert to list class
#' @export
#' @param x Object of class hdpState
#' @param ... unused
setMethod("as.list",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- list(numdp=as.numeric(x@numdp),
                        numconparam=x@numconparam,
                        base=as.list(x@base),
                        conparam=lapply(x@conparam, as.list),
                        dp=lapply(x@dp, as.list),
                        dpstate=x@dpstate,
                        ppindex=x@ppindex,
                        cpindex=x@cpindex,
                        ttindex=x@ttindex)
            return(ans)
          })

# convert list to hdpState class
setGeneric("as.hdpState", function(x, ...) standardGeneric("as.hdpState"))
setMethod("as.hdpState",
          signature = "list",
          definition = function(x, ...) {
            ans <- new("hdpState",
                       numdp=as.integer(x$numdp),
                       numconparam=x$numconparam,
                       base=as.hdpBase(x$base),
                       conparam=lapply(x$conparam, as.hdpConparam),
                       dp=lapply(x$dp, as.hdpDP),
                       dpstate=x$dpstate,
                       ppindex=x$ppindex,
                       cpindex=x$cpindex,
                       ttindex=x$ttindex)
            return(ans)
          })


# show method
setMethod("show",
          "hdpState",
          function(object) {
            cat("Object of class", class(object), "\n")
            numdp <- object@numdp
            cat(" Number of DP nodes:", numdp, "\n")
            cat(" Index of parent DP:",
                object@ppindex[1:min(10, numdp)], "...\n")
            cat(" Number of data items per DP:",
                sapply(object@dp[1:min(10, numdp)],
                       function(x) x@numdata), "...\n")
            cat(" Index of conparam per DP:",
                object@cpindex[1:min(10, numdp)], "...\n")
            cat(" Conparam hyperparameters and current value:\n")
            print(matrix(sapply(object@conparam,
                                function(x) c(x@alphaa, x@alphab, x@alpha)),
                         nrow=length(object@conparam),
                         byrow=T,
                         dimnames=list(
                           paste("Conparam", 1:length(object@conparam)),
                           c("Shape", "Rate", "Value"))))
            cat(" Number of data categories:", length(object@base@hh), "\n")
            cat(" Number of clusters:", object@base@numclass, "\n")
          })


# getter methods for hdpState objects ---------

setGeneric("numdp", function(x, ...) standardGeneric("numdp"))
#' @describeIn hdpState Get number of DPs
#' @aliases numdp
#' @export
setMethod("numdp",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@numdp
            return(ans)
          })

setGeneric("numconparam", function(x, ...) standardGeneric("numconparam"))
#' @describeIn hdpState Get number of concentration parameters
#' @aliases numconparam
#' @export
setMethod("numconparam",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@numconparam
            return(ans)
          })

setGeneric("base", function(x, ...) standardGeneric("base"))
#' @describeIn hdpState Get base distribution
#' @aliases base
#' @export
setMethod("base",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@base
            return(ans)
          })

setGeneric("conparam", function(x, ...) standardGeneric("conparam"))
#' @describeIn hdpState Get list of concentration parameters
#' @aliases conparam
#' @export
setMethod("conparam",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@conparam
            return(ans)
          })

setGeneric("dp", function(x, ...) standardGeneric("dp"))
#' @describeIn hdpState Get list of DP nodes
#' @aliases dp
#' @export
setMethod("dp",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@dp
            return(ans)
          })

setGeneric("dpstate", function(x, ...) standardGeneric("dpstate"))
#' @describeIn hdpState Get state of every DP
#' @aliases dpstate
#' @export
setMethod("dpstate",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@dpstate
            return(ans)
          })

setGeneric("ppindex", function(x, ...) standardGeneric("ppindex"))
#' @describeIn hdpState Get parent process index of every DP
#' @aliases ppindex
#' @export
setMethod("ppindex",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@ppindex
            return(ans)
          })

setGeneric("cpindex", function(x, ...) standardGeneric("cpindex"))
#' @describeIn hdpState Get concentration parameter index of every DP
#' @aliases cpindex
#' @export
setMethod("cpindex",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@cpindex
            return(ans)
          })

setGeneric("numcateg", function(x, ...) standardGeneric("numcateg"))
#' @describeIn hdpState Get number of data categories
#' @aliases numcateg
#' @export
setMethod("numcateg",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- length(x@base@hh)
            return(ans)
          })


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
                        clust_dp_weights=x@clust_dp_weights,
                        comp_settings=x@comp_settings,
                        comp_categ_counts=x@comp_categ_counts,
                        comp_dp_counts=x@comp_dp_counts,
                        comp_dp_weights=x@comp_dp_weights,
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
            cat(" Conparam moves per iter:", object@settings$cpiter, "\n")
            cat(" Tabulate raw cluster number across samples:")
            print(table(object@numcluster))
            ncomp <- length(object@comp_categ_counts)
            if (ncomp == 0) {
              cat(" Components: NO. Run hdp_extract_components \n")
            } else {
              cat(" Components: YES. Prop of data explained >",
                  object@comp_settings$prop.ex,
                  " Merge if cosine sim >",
                  object@comp_settings$cos.merge, "\n")
              cat(" Component number:", ncomp-1, "\n")
            }
            cat(" ----------\n")
            cat(" Final hdpState: \n")
            print(object@hdp)
          })


# getter methods for hdpSampleChain objects ---------

setGeneric("hdp_seed", function(x, ...) standardGeneric("hdp_seed"))
#' @describeIn hdpSampleChain Get random seed of posterior sampling chain
#' @aliases hdp_seed
#' @export
setMethod("hdp_seed",
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

setGeneric("clust_dp_weights",
           function(x, ...) standardGeneric("clust_dp_weights"))
#' @describeIn hdpSampleChain Get dp node vs cluster weights for each posterior sample
#' @aliases clust_dp_weights
#' @export
setMethod("clust_dp_weights",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@clust_dp_weights
            return(ans)
          })

setGeneric("hdp_comp_settings",
           function(x, ...) standardGeneric("hdp_comp_settings"))
#' @describeIn hdpSampleChain Get settings of component extraction from HDP posterior sampling chain
#' @aliases hdp_comp_settings
#' @export
setMethod("hdp_comp_settings",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@comp_settings
            return(ans)
          })


setGeneric("comp_categ_counts",
           function(x, ...) standardGeneric("comp_categ_counts"))
#' @describeIn hdpSampleChain Get sample vs category counts for each component
#' @aliases comp_categ_counts
#' @export
setMethod("comp_categ_counts",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@comp_categ_counts
            return(ans)
          })

setGeneric("comp_dp_counts",
           function(x, ...) standardGeneric("comp_dp_counts"))
#' @describeIn hdpSampleChain Get sample vs component counts for each DP
#' @aliases comp_dp_counts
#' @export
setMethod("comp_dp_counts",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@comp_dp_counts
            return(ans)
          })

setGeneric("comp_dp_weights",
           function(x, ...) standardGeneric("comp_dp_weights"))
#' @describeIn hdpSampleChain Get sample vs component weights for each DP
#' @aliases comp_dp_weights
#' @export
setMethod("comp_dp_weights",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@comp_dp_weights
            return(ans)
          })

setGeneric("comp_categ_distn",
           function(x, ...) standardGeneric("comp_categ_distn"))
#' @describeIn hdpSampleChain Get mean distribution over data categories for each component
#' @aliases comp_categ_distn
#' @export
setMethod("comp_categ_distn",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@comp_categ_distn
            return(ans)
          })

setGeneric("comp_dp_distn",
           function(x, ...) standardGeneric("comp_dp_distn"))
#' @describeIn hdpSampleChain Get mean distribution over components for each DP
#' @aliases comp_dp_distn
#' @export
setMethod("comp_dp_distn",
          signature = "hdpSampleChain",
          definition = function(x, ...) {
            ans <- x@comp_dp_distn
            return(ans)
          })
