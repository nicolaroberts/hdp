#' Diagnostic plots for HDP posterior sampling chain
#'
#' @param chain A hdpSampleChain object
#' @name plotchain
#' @examples
#' plot_lik(mut_example_chain, bty="L")
#' plot_numcluster(mut_example_chain, bty="L")
#' plot_data_assigned(mut_example_chain, bty="L")
NULL
#> NULL

#' @param start The starting iteration to plot from (default 1)
#' @param end The final iteration to plot to (default is end of chain)
#' @param col_lik Plot colour of likelihood (default blue)
#' @param col_burn Plot colour of burnin (default red)
#' @param xlab Horizontal axis label
#' @param ylab Vertical axis label
#' @param ... Other arguments to plot
#' @export
#' @rdname plotchain
plot_lik <- function(chain, start=1, end=length(lik(chain)),
                     col_lik="blue", col_burn="red",
                     xlab="Iteration", ylab="Likelihood", ...){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (length(start) != 1 | start < 1 | start > length(lik(chain)) |
        start %% 1 != 0 | length(end) != 1 | end < 1 |
        end > length(lik(chain)) | end %% 1 != 0 ) {
    stop("start and end must be positive integers no greater than iter number")
  }

  # plot - max 5000 points
  iter <- ceiling(seq(start, end, length=min(end - start, 5000)))
  lik <- lik(chain)[iter]
  plot(iter, lik, type="l", col=col_lik,
       xlab=xlab, ylab=ylab, ...)
  abline(v=hdp_settings(chain)$burnin, col=col_burn, lty=2)
}

#' @param col Plot colour for numcluster (default blue)
#' @export
#' @rdname plotchain
plot_numcluster <- function(chain, col="blue", xlab="Sample",
                            ylab="Number of raw clusters", ...){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")

  # plot
  numcluster <- numcluster(chain)
  plot(numcluster, type="l", col=col, xlab=xlab, ylab=ylab, ...)
}


#' @param legend Logical - should a legend be included? (default TRUE)
#' @param col_early Color ramp side for early posterior samples
#' @param col_late Color ramp side for late posterior samples
#' @param dat_prop Extend horiztonal axis to dat_prop proportion of data assigned
#' @export
#' @rdname plotchain
plot_data_assigned <- function(chain, legend=TRUE, col_early="hotpink",
                               col_late="skyblue3",
                               dat_prop = 0.995,
                               xlab="Number of raw clusters",
                               ylab="Cumulative prop. of data assigned", ...){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if(class(legend) != "logical") stop("legend must be TRUE or FALSE")


  # extract necessary info
  n <- hdp_settings(chain)$n
  ccc <- clust_categ_counts(chain)
  numcluster <- numcluster(chain)

  # func to adjust length of a vector
  set.length <- function(vec, len){
    length(vec) <- len
    return(vec)
  }

  # calc cumulative proportion of data items assigned to a cluster (per sample),
  # ordering each cluster from largest to smallest
  pda <- sapply(ccc, function(x) cumsum(sort(colSums(x), decreasing=T))/sum(x))
  pda <- sapply(pda, set.length, len=max(numcluster))

  # average line
  avg <- rowMeans(pda, na.rm=T)
  # number of clusters (on average) needed to explain 99.9% of dataset
  xmax <- which(avg>dat_prop)[1]

  # colour ramp across posterior samples
  cols <- colorRampPalette(colors=c(col_early, col_late))

  matplot(pda[1:xmax,], type="l", col=cols(n), xlim=c(0.95, xmax),
          xlab=xlab, ylab=ylab, ...)

  if (legend) {
    legend("bottomright", col=c(col_early, col_late, "black"), lty=1, lwd=2,
           legend=c("Early samples", "Late samples", "Average"), bty="n")
  }

  points(avg[1:xmax], type="l", cex=2)

}
