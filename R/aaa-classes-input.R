#' hdpBase class for the base distribution
#' @export
#'
#' @slot hh parameters for base Dirichlet distribution (pseudocounts)
#' @slot classqq overall count matrix for data items of each category (rows) in each cluster (columns)
#' @slot numclass number of clusters
setClass("hdpBase",
         slots = list(
           hh = "numeric",
           classqq = "matrix",
           numclass = "integer"),
         validity = function(object) {
           is_valid <- TRUE
           if (is.any.slot.negative(object)) {
             is_valid <- FALSE
             message("hdpBase must not contain negative values")
           }
           if (nrow(object@classqq) != length(object@hh)) {
             is_valid <- FALSE
             message("nrow(classqq) must equal length(hh)")
           }
           if (ncol(object@classqq) != (object@numclass + 1)) {
             is_valid <-  FALSE
             message("ncol(classqq) must equal numclass+1")
           }
           return(is_valid)
         })


#' hdpConparam class for the DP concentration parameter/s
#' @export
#'
#' @slot alphaa shape parameter for the gamma prior over alpha
#' @slot alphab rate parameter for the gamma prior over alpha
#' @slot numdp number of DP nodes sharing this concentration parameter
#' @slot alpha concentration parameter value
#' @slot totalnd number of data items in each DP with this concentration parameter
#' @slot totalnt number of tables in each DP with this concentration parameter
setClass("hdpConparam",
         slots = list(
           alphaa = "numeric",
           alphab = "numeric",
           numdp = "integer",
           alpha = "numeric",
           totalnd = "integer",
           totalnt = "integer"),
         validity = function(object){
           is_valid <- TRUE
           if (is.any.slot.negative(object)) {
             is_valid <- FALSE
             message("hdpConparam must not contain negative values")
           }
           if (length(object@totalnd) != object@numdp |
                 length(object@totalnt) != object@numdp) {
             is_valid <- FALSE
             message("length(totalnd/nt) must be numdp in conparam object")
           }
           return(is_valid)
         })


#' hdpDP class for a DP node
#'
#' note that the 'items' in parent nodes are the tables of their children
#'
#' @export
#'
#' @slot datacc cluster index for each data item
#' @slot classnd number of items assigned to each cluster in this DP
#' @slot classnt number of tables assigned to each cluster in this DP
#' @slot beta weight on each cluster in this DP (including empty cluster at end)
#' @slot alpha concentration parameter for this DP
#' @slot numdata number of data items registered to this DP node
#' @slot datass value of each data item
setClass("hdpDP",
         slots = list(
           datacc = "integer",
           classnd = "integer",
           classnt = "integer",
           beta = "numeric",
           alpha = "numeric",
           numdata = "integer",
           datass = "integer"),
         validity = function(object){
           is_valid <- TRUE
           if (is.any.slot.negative(object)) {
             is_valid <- FALSE
             message("hdpDP must not contain negative values")
           }
           if (length(object@datass) != object@numdata |
                 length(object@datacc) > 0 &
                 length(object@datacc) != object@numdata) {
             is_valid <- FALSE
             message("length of datass and datacc must be numdata")
           }
           return(is_valid)
         })



#' hdpState class for a Hierarchical Dirichlet Process in one state
#'
#' @slot numdp number of DP nodes in the hierarchical Dirichlet Process
#' @slot numconparam number of concentration parameters
#' @slot base base distribution (hdpBase object)
#' @slot conparam concentration parameters (list of hdpConparam objects)
#' @slot dp DP nodes (list of hdpDP objects)
#' @slot dpstate state of DP nodes for posterior sampling process: active (2), frozen (1), or heldout (0)
#' @slot ppindex parent node index for each DP
#' @slot cpindex concentration parameter index for each DP
#' @slot ttindex DP index of those sharing a concentration parameter
#' @slot initcc number of initial clusters
#' @slot seed_activate random seed used to initiate cluster membership
#' @slot pseudoDP (Optional) index of pseudodata nodes (only if initialised via hdp_prior_init)
#' @export
setClass("hdpState",
         slots = list(
           numdp = "integer",
           numconparam = "integer",
           base = "hdpBase",
           conparam = "list",
           dp = "list",
           dpstate = "integer",
           ppindex = "integer",
           cpindex = "integer",
           ttindex = "integer",
           initcc = "integer",
           seed_activate = "integer",
           pseudoDP = "integer"),
         validity = function(object){
           is_valid <- TRUE
           if (is.any.slot.negative(object)) {
             is_valid <- FALSE
             message("hdpState must not contain negative values")
           }
           if (!validObject(object@base)) {
             is_valid <- FALSE
             message("base not valid")
           }
           if (object@numconparam > 0 &
                 any(!sapply(object@conparam, validObject))) {
             is_valid <- FALSE
             message("conparam not valid")
           }
           if (object@numdp > 0 & any(!sapply(object@dp, validObject))) {
             is_valid <- FALSE
             message("dp not valid")
           }
           if (length(object@conparam) != object@numconparam) {
             is_valid <- FALSE
             message("length(conparam) must be numconparam")
           }
           if (length(object@dp) != object@numdp) {
             is_valid <- FALSE
             message("length(dp) must be numdp")
           }
           if (!any(object@dpstate %in% c(0,1,2))) {
             is_valid <- FALSE
             message("dpstate must only have values 0, 1, 2")
           }
           if (any(object@ppindex > object@numdp) |
                 any(object@ttindex > object@numdp)) {
             is_valid <- FALSE
             message("ppindex and ttindex can not be greater than numdp")
           }
           if(any(object@cpindex > object@numconparam)) {
             is_valid <- FALSE
             message("cpindex can not be greater than numconparam")
           }
           return(is_valid)
         })
