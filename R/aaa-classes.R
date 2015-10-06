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
#' @slot numdp number of DPs sharing this concentration parameter
#' @slot alpha concentration parameter value
#' @slot totalnd number of items in each DP with this concentration parameter
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
#' @slot numconparam number of different concentration parameters
#' @slot base base distribution (hdpBase object)
#' @slot conparam concentration parameters (list of hdpConparam objects)
#' @slot dp DP nodes (list of hdpDP objects)
#' @slot dpstate state of DP nodes for posterior sampling process: active (2), frozen (1), or heldout (0)
#' @slot ppindex parent node index for each DP
#' @slot cpindex concentration parameter index for each DP
#' @slot ttindex DP index of those sharing a concentration parameter
#' @slot initcc Number of initial clusters
#' @slot seed_activate Random seed used to initiate cluster membership
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
           seed_activate = "integer"),
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


#' hdpSampleChain class for posterior samples of a hierarchical Dirichlet Process
#'
#' @export
#'
#' @slot seed Random seed
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
#' @slot clust_dp_weights List of matrices (one from each posterior sample)
#'  with the weights of each cluster per DP. Number of rows is the number of
#'  DPs (constant), and number of columns is the number of clusters in that
#'  posterior sample (variable).
#' @slot comp_settings Settings used to
#'  consolidate raw clusters (variable number across posterior samples) into
#'  global components (constant number across posterior samples):
#'    prop.ex (minimum proportion of data explained), cos.merge (merge components within cosine similarity)
#' @slot comp_categ_counts List of matrices (one for each component)
#'  counting the sample-category data assignment across all DP nodes.
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of data categories.
#' @slot comp_dp_counts List of matrices (one for each DP)
#'  counting sample-component assignment (aggregating across data categories).
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of components.
#' @slot comp_dp_weights List of matrices (one for each DP)
#'  with the weights of each component per sample. Number of rows is the number of
#'  posterior samples, and number of columns is the number of components.
#' @slot comp_categ_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95% credibility interval) distribution
#'  over data categories for each component. Number of rows is the number of
#'  components, and number of columns is the number of data categories. Rows sum to 1.
#' @slot comp_dp_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95% credibility interval) distribution
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
           clust_dp_weights = "list",
           comp_settings = "list",
           comp_categ_counts = "list",
           comp_dp_counts = "list",
           comp_dp_weights="list",
           comp_categ_distn="list",
           comp_dp_distn="list"),
         validity = function(object){
           is_valid <- TRUE
           if (!validObject(object@hdp)){
             is_valid <- FALSE
             message("final hdpState not valid")
           }
           nsamp <- length(object@numcluster)
           if (nrow(object@cp_values) != nsamp |
                 length(object@clust_categ_counts) != nsamp |
                 length(object@clust_dp_counts) != nsamp |
                 length(object@clust_dp_weights) != nsamp) {
             is_valid <- FALSE
             message("inconsistent sample number")
           }
           nclust <- object@numcluster + 1
           if (any(sapply(object@clust_categ_counts, ncol) != nclust) |
                 any(sapply(object@clust_dp_counts, ncol) != nclust) |
                 any(sapply(object@clust_dp_weights, ncol) != nclust)) {
             is_valid <- FALSE
             message("inconsistent cluster number")
           }
           ncateg <- length(object@hdp@base@hh)
           if (any(sapply(object@clust_categ_counts, nrow) != ncateg)) {
             is_valid <- FALSE
             message("inconsistent category number")
           }
           if (any(signif(sapply(object@clust_dp_weights, rowSums), 3) != 1)){
             is_valid <- FALSE
             message("clust_dp_weights do not all sum to 1")
           }
           settings <- sapply(object@settings, function(x) x)
           if (!all.equal(names(settings), c("burnin","n","space","cpiter"))){
             is_valid <- FALSE
             message("settings list incorrectly named")
           }
           if (any(settings %% 1 != 0) | any(settings < 1)){
             is_valid <- FALSE
             message("settings must be positive integers")
           }
           # add validity checks for comp* slots and seeds list
           return(is_valid)
         })



#' hdpSampleMulti class for multiple independent hdpSampleChain objects for the same HDP
#'
#' @export
#'
#' @slot chains List of hdpSampleChain objects storing multiple independent runs of the posterior sampling chain for the same data and HDP struct
#' @slot comp_settings Settings used to
#'  consolidate raw clusters (variable number across posterior samples) into
#'  global components (constant number across posterior samples):
#'    cos.merge (merge components within cosine similarity)
#' @slot comp_categ_counts List of matrices (one for each component)
#'  counting the sample-category data assignment across all DP nodes.
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of data categories.
#' @slot comp_dp_counts List of matrices (one for each DP)
#'  counting sample-component assignment (aggregating across data categories).
#'  Number of rows is the number of posterior samples, and number of
#'  columns is the number of components.
#' @slot comp_dp_weights List of matrices (one for each DP)
#'  with the weights of each component per sample. Number of rows is the number of
#'  posterior samples, and number of columns is the number of components.
#' @slot comp_categ_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95% credibility interval) distribution
#'  over data categories for each component. Number of rows is the number of
#'  components, and number of columns is the number of data categories. Rows sum to 1.
#' @slot comp_dp_distn List with elements "mean" and "cred.int", containing
#'  matrices with the mean (and lower/upper 95% credibility interval) distribution
#'  over components for each DP. Number of rows is the number of
#'  DPs, and number of columns is the number of components. Rows sum to 1.
#' @slot prop.ex (Average) proportion of dataset explained by the components consistently
#'  found in all chains.
#'
setClass("hdpSampleMulti",
         slots = list(
           chains = "list",
           comp_settings = "list",
           comp_categ_counts = "list",
           comp_dp_counts = "list",
           comp_dp_weights="list",
           comp_categ_distn="list",
           comp_dp_distn="list",
           prop.ex="numeric"),
         validity = function(object){
           is_valid <- TRUE
           # add validity checks
           # check all hdpSampleChain objects have same data, struct
           # and different seed.
           return(is_valid)
         })
