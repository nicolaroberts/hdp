#' hdpSampleChain class for posterior samples off one MCMC chain
#'
#' @export
#'
#' @slot seed Random seed used by \code{hdp_posterior}
#' @slot settings Settings of the posterior sampling chain: burnin, n (number of samples collected),
#'  space (iters between samples), cpiter (con param moves between iters)
#' @slot hdp hdpState object after the final iteration
#' @slot lik Likelihood of data given model at each iteration
#' @slot numcluster Number of raw data clusters in each posterior sample
#' @slot cp_values Matrix of concentration parameter values (one column for each parameter) in each posterior sample (rows).
#' @slot clust_categ_counts List of matrices (one from each posterior sample)
#'  counting the category-cluster data assignment across all DP nodes.
#'  Number of rows is the number of categories (constant), and number of
#'  columns is the number of clusters in that posterior sample (variable).
#' @slot clust_dp_counts List of matrices (one from each posterior sample)
#'  counting within-DP cluster assignment (aggregating across data categories).
#'  Number of rows is the number of DPs (constant), and number of
#'  columns is the number of clusters in that posterior sample (variable).
#' @slot numcomp Number of global components extracted by \code{hdp_extract_components}
#'  (not including component 0)
#' @slot prop.ex (Average) proportion of dataset explained by the extracted components
#' @slot comp_cos_merge \code{cos.merge} setting used by \code{hdp_extract_components}
#' @slot comp_categ_counts List of matrices (one for each component)
#'  counting the sample-category data assignment across all DP nodes.
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of data categories.
#' @slot comp_dp_counts List of matrices (one for each DP)
#'  counting sample-component assignment (aggregating across data categories).
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of components.
#' @slot comp_categ_distn List with elements \code{mean} and \code{cred.int}, containing
#'  matrices with the mean (and lower/upper 95\% credibility interval) distribution
#'  over data categories for each component. Number of rows is the number of
#'  components, and number of columns is the number of data categories. Rows sum to 1.
#' @slot comp_dp_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95\% credibility interval) distribution
#'  over components for each DP. Number of rows is the number of
#'  DPs, and number of columns is the number of components. Rows sum to 1.
#'
setClass("hdpSampleChain",
         slots = list(
           seed = "integer",
           settings = "list",
           hdp = "hdpState",
           lik = "numeric",
           numcluster = "integer",
           cp_values = "matrix",
           clust_categ_counts = "list",
           clust_dp_counts = "list",
           numcomp = "integer",
           prop.ex = "numeric",
           comp_cos_merge = "numeric",
           comp_categ_counts = "list",
           comp_dp_counts = "list",
           comp_categ_distn="list",
           comp_dp_distn="list"),
         validity = function(object){
           is_valid <- TRUE

           if (!validObject(object@hdp)){
             is_valid <- FALSE
             message("final hdpState not valid")
           }

           settings <- sapply(object@settings, function(x) x)
           if (!all.equal(names(settings), c("burnin","n","space","cpiter"))){
             is_valid <- FALSE
             message("settings list incorrectly named")
           }
           if (any(settings %% 1 != 0) || any(settings < 1)){
             is_valid <- FALSE
             message("settings must be positive integers")
           }

           is_comp <- length(object@numcomp)==1

           nsamp <- object@settings$n
           if (length(object@numcluster) != nsamp ||
               nrow(object@cp_values) != nsamp ||
               length(object@clust_categ_counts) != nsamp ||
               length(object@clust_dp_counts) != nsamp ||
               is_comp && any(sapply(object@comp_categ_counts, nrow) != nsamp) ||
               is_comp && any(sapply(object@comp_dp_counts, nrow) != nsamp)) {
             is_valid <- FALSE
             message("inconsistent sample number")
           }

           nclust <- object@numcluster + 1
           if (any(sapply(object@clust_categ_counts, ncol) != nclust) ||
               any(sapply(object@clust_dp_counts, ncol) != nclust)) {
             is_valid <- FALSE
             message("inconsistent cluster number")
           }

           ncateg <- length(object@hdp@base@hh)
           if (any(sapply(object@clust_categ_counts, nrow) != ncateg) ||
               is_comp && any(sapply(object@comp_categ_counts, ncol) != ncateg) ||
               is_comp && ncol(object@comp_categ_distn$mean) != ncateg) {
             is_valid <- FALSE
             message("inconsistent category number")
           }

           ndp <- object@hdp@numdp
           if (any(sapply(object@clust_dp_counts, nrow) != ndp) ||
               is_comp && length(object@comp_dp_counts) != ndp ||
               is_comp && nrow(object@comp_dp_distn$mean) != ndp) {
             is_valid <- FALSE
             message("inconsistent DP number")
           }

           if (is_comp){
             if (object@prop.ex <= 0 ||
                 object@prop.ex > 1) {
               is_valid <- FALSE
               message("prop.ex must be between 0 and 1")
             }

             if (object@comp_cos_merge <= 0 ||
                 object@comp_cos_merge > 1) {
               is_valid <- FALSE
               message("comp_cos_merge must be between 0 and 1")
             }

             ncomp <- object@numcomp + 1
             if (length(object@comp_categ_counts) != ncomp ||
                 any(sapply(object@comp_dp_counts, ncol) != ncomp) ||
                 nrow(object@comp_categ_distn$mean) != ncomp ||
                 ncol(object@comp_dp_distn$mean) != ncomp) {
               is_valid <- FALSE
               message("inconsistent component number")
             }
           }
           return(is_valid)
         })



#' hdpSampleMulti class for multiple independent hdpSampleChain objects for the same HDP
#'
#' @export
#'
#' @slot chains List of hdpSampleChain objects storing multiple independent runs of the posterior sampling chain for the same data and HDP struct
#' @slot numcomp Number of global components extracted by \code{hdp_extract_components}
#'  (not including component 0)
#' @slot prop.ex (Average) proportion of dataset explained by the extracted components
#' @slot comp_cos_merge \code{cos.merge} setting used by \code{hdp_extract_components}
#' @slot comp_categ_counts List of matrices (one for each component)
#'  counting the sample-category data assignment across all DP nodes.
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of data categories.
#' @slot comp_dp_counts List of matrices (one for each DP)
#'  counting sample-component assignment (aggregating across data categories).
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of components.
#' @slot comp_categ_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95\% credibility interval) distribution
#'  over data categories for each component. Number of rows is the number of
#'  components, and number of columns is the number of data categories. Rows sum to 1.
#' @slot comp_dp_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95\% credibility interval) distribution
#'  over components for each DP. Number of rows is the number of
#'  DPs, and number of columns is the number of components. Rows sum to 1.
#'
setClass("hdpSampleMulti",
         slots = list(
           chains = "list",
           numcomp = "integer",
           prop.ex="numeric",
           comp_cos_merge = "numeric",
           comp_categ_counts = "list",
           comp_dp_counts = "list",
           comp_categ_distn="list",
           comp_dp_distn="list"),
         validity = function(object){
           is_valid <- TRUE

           # check all constituent hdpSampleChain objects are valid,
           # with same data and struct, and different seeds
           if (any(!sapply(object@chains, validObject))){
             is_valid <- FALSE
             message("constituent hdpSampleChain not valid")
           }

           all_hdpStates <- sapply(object@chains, final_hdpState)
           if (any(diff(sapply(all_hdpStates, slot, "numdp")) != 0) ||
               any(diff(sapply(all_hdpStates, slot, "numconparam")) != 0) ||
               any(apply(sapply(all_hdpStates, slot, "dpstate"), 1, diff) != 0) ||
               any(apply(sapply(all_hdpStates, slot, "ppindex"), 1, diff) != 0) ||
               any(apply(sapply(all_hdpStates, slot, "cpindex"), 1, diff) != 0) ||
               any(diff(sapply(all_hdpStates, function(x) sum(x@base@classqq))) != 0) ||
               any(diff(sapply(all_hdpStates, function(x) length(x@base@hh))) != 0)) {
             is_valid <- FALSE
             message("constituent hdpState structures not identical")
           }

           if (length(unique(sapply(object@chains, slot, "seed"))) != length(object@chains)) {
             is_valid <- FALSE
             message("constituent hdpSampleChains must have different random seeds")
           }

           is_comp <- length(object@numcomp)==1

           if (is_comp) {

             if (object@prop.ex <= 0 ||
                 object@prop.ex > 1) {
               is_valid <- FALSE
               message("prop.ex must be between 0 and 1")
             }

             if (object@comp_cos_merge <= 0 ||
                 object@comp_cos_merge > 1) {
               is_valid <- FALSE
               message("comp_cos_merge must be between 0 and 1")
             }

             nsamp <- sum(sapply(object@chains, function(x) x@settings$n))
             if (any(sapply(object@comp_categ_counts, nrow) != nsamp) ||
                 any(sapply(object@comp_dp_counts, nrow) != nsamp)) {
               is_valid <- FALSE
               message("inconsistent sample number")
             }

             ndp <- object@chains[[1]]@hdp@numdp
             if (length(object@comp_dp_counts) != ndp ||
                 nrow(object@comp_dp_distn$mean) != ndp) {
               is_valid <- FALSE
               message("inconsistent DP number")
             }

             ncomp <- object@numcomp + 1
             if (length(object@comp_categ_counts) != ncomp ||
                 any(sapply(object@comp_dp_counts, ncol) != ncomp) ||
                 nrow(object@comp_categ_distn$mean) != ncomp ||
                 ncol(object@comp_dp_distn$mean) != ncomp) {
               is_valid <- FALSE
               message("inconsistent component number")
             }
           }

           return(is_valid)
         })
