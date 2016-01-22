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
            cat(" Prior over data categories:", object@hh[1:min(5, numcat)], "...\n")
            cat(" Number of raw clusters:", numclass, "\n")
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
            if (object@numdata == 0){
              cat(" Number of child data tables in each cluster:",
                  object@classnd, "\n")
            } else {
              cat(" Number of data items in each cluster:",
                  object@classnd, "\n")
              cat(" Number of data items in each category:",
                  table(object@datass), "\n")
            }
          })

setGeneric("numdata", function(x, ...) standardGeneric("numdata"))
#' @describeIn hdpDP Get number of data items at this DP.
#' @aliases numdata
#' @export
setMethod("numdata",
          signature = "hdpDP",
          definition = function(x, ...) {
            ans <- x@numdata
            return(ans)
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
                        ttindex=x@ttindex,
                        initcc=x@initcc,
                        seed_activate=x@seed_activate,
                        pseudoDP=x@pseudoDP)
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
                       ttindex=x$ttindex,
                       initcc=x$initcc,
                       seed_activate=x$seed_activate,
                       pseudoDP=x$pseudoDP)
            return(ans)
          })

# show method
setMethod("show",
          "hdpState",
          function(object) {
            cat("Object of class", class(object), "\n")
            numdp <- object@numdp
            cat(" Number of DP nodes:", numdp, "\n")

            npseudo <- length(object@pseudoDP)
            if (npseudo > 0) {
              cat(" Number of pseudo DP nodes for prior info:",
                  npseudo, "\n")
              cat(" Index of pseudo DP nodes:",
                  object@pseudoDP[1:min(10, npseudo)], "...\n")
            }
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
            if (length(object@initcc) == 1) {
              cat(" Initialised with", object@initcc,
                  "clusters, using random seed", object@seed_activate, "\n")
            }
          })

# getter methods for hdpState objects

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

setGeneric("base_params", function(x, ...) standardGeneric("base_params"))
#' @describeIn hdpState Get parameters of the base Dirichlet distribution (like psuedocounts across categories).
#' @aliases base_params
#' @export
setMethod("base_params",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@base@hh
            return(ans)
          })

setGeneric("activating_seed", function(x, ...) standardGeneric("activating_seed"))
#' @describeIn hdpState Get seed used to initialse clustering
#' @aliases activating_seed
#' @export
setMethod("activating_seed",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@seed_activate
            return(ans)
          })

setGeneric("pseudoDP", function(x, ...) standardGeneric("pseudoDP"))
#' @describeIn hdpState Get index of frozen pseudo-data DP nodes for prior info
#'  (only if initialised via hdp_prior_init)
#' @aliases pseudoDP
#' @export
setMethod("pseudoDP",
          signature = "hdpState",
          definition = function(x, ...) {
            ans <- x@pseudoDP
            return(ans)
          })
