#' Assign data to DP nodes in a hdpState object
#'
#' Assign data to 'heldout' (state is 0) DP nodes in a hdpState object.
#' 'Heldout' DPs are not available for posterior sampling,
#' and will need to be activated (see \code{\link{dp_activate}}). The posterior
#' sampling process (a Gibbs sampler) is run via \code{\link{hdp_posterior}}.
#'
#' @param hdp A hdpState object
#' @param dpindex Indices of the DPs to assign data to (in same order as rows of \code{data})
#' @param data A \code{data.frame} or \code{matrix} of counts with one row for every sample
#' (same order as \code{dpindex}) and one column for every data category.
#' @return A hdpState object updated with the new data. See \code{\link{hdpState-class}}
#' @seealso \code{\link{hdp_init}}, \code{\link{hdp_adddp}},
#'  \code{\link{dp_activate}}, \code{\link{hdp_posterior}}
#' @export
#' @examples
#' example_data_hdp
#' my_hdp <- hdp_init(ppindex=0, cpindex=1, hh=rep(1, 6), alphaa=rep(1, 3), alphab=rep(2, 3))
#' my_hdp <- hdp_adddp(my_hdp, 2, 1, 2)
#' my_hdp <- hdp_adddp(my_hdp, 10, c(rep(2, 5), rep(3, 5)), 3)
#' my_hdp <- hdp_setdata(my_hdp, 4:13, example_data_hdp)
#' dp(my_hdp)

hdp_setdata <- function(hdp, dpindex, data){

  #input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  if (any(dpindex < 1) |
        any(dpindex > hdp@numdp) |
        any(dpindex %% 1 != 0) |
        any(duplicated(dpindex))) {
    stop("dpindex must be positive integers no greater than
         numdp(hdp) with no duplicates")
  }
  if (!class(data) %in% c("matrix", "data.frame")) {
    stop("data must be data.frame or matrix")
  }
  if (nrow(data)!=length(dpindex)) stop("nrow(data) must equal length(dpindex)")
  if (ncol(data)!=numcateg(hdp)) stop("ncol(data) must equal numcateg(hdp)")
  if (any(data %% 1 != 0) | any(data < 0)) {
    stop("data must contain non-negative integer values")
  }

  # convert data to a list of data item values (not category counts)
  datass <- apply(data, 1, function(x) rep(1:ncol(data), x))


  HELDOUT <- 0L

  # assign data to specified DP
  for (jj in 1:length(dpindex)){
    if (hdp@dpstate[dpindex[jj]] != HELDOUT){
      stop("Cannot set data for DPs that are not held out")
    }
    hdp@dp[[dpindex[jj]]]@numdata <- length(datass[[jj]])
    hdp@dp[[dpindex[jj]]]@datass  <- datass[[jj]]
  }

  # check validity and return
  if (!validObject(hdp)) warning("Not a valid hdpState object.")
  return(hdp)
}
