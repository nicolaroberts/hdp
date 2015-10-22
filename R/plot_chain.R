#' Plot likelihood of HDP posterior sampling chain.
#'
#' @param chain A hdpSampleChain object
#' @param start The starting iteration to plot from (default 1)
#' @param end The final iteration to plot to (default is end of chain)
#' @param col_lik Plot colour of likelihood (default blue)
#' @param col_burn Plot colour of burnin (default red)
#' @param xlab Horizontal axis label
#' @param ylab Vertical axis label
#' @param ... Other arguments to plot
#' @seealso \code{\link{hdpSampleChain-class}}, \code{\link{plot_numcluster}},
#'  \code{\link{plot_data_assigned}}
#' @export
#' @examples
#' plot_lik(mut_example_chain, bty="L")
#'

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


#' Plot number of raw clusters in each posterior sample from the chain
#'
#' @param chain A hdpSampleChain object
#' @param col Plot colour (default blue)
#' @param xlab Horizontal axis label
#' @param ylab Vertical axis label
#' @param ... Other arguments to plot
#' @seealso \code{\link{hdpSampleChain-class}}, \code{\link{plot_lik}},
#'  \code{\link{plot_data_assigned}}
#' @export
#' @examples
#' plot_numcluster(mut_example_chain, bty="L")
#'
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


#' Plot proportion of data items assigned vs number of clusters
#'
#' Plot the cumulative proportion of data items assigned as the number of raw
#' clusters increases (raw clusters sorted from largest to smallest).
#' The horizontal axis extends as far as needed to explain 99.9\% of data
#' items (on average).
#' Each posterior sample is plotted separately, using a colour spectrum ranging from
#' \code{col_early} to \code{col_late}. The average across posterior samples
#' is plotted in black.
#'
#' @param chain A hdpSampleChain object
#' @param legend Logical - should a legend be included? (default TRUE)
#' @param col_early Color ramp side for early posterior samples
#' @param col_late Color ramp side for late posterior samples
#' @param xlab Horizontal axis label
#' @param ylab Vertical axis label
#' @param ... Other arguments to plot
#' @seealso \code{\link{hdpSampleChain-class}}, \code{\link{plot_lik}},
#'  \code{\link{plot_numcluster}}
#' @export
#' @examples
#' plot_data_assigned(mut_example_chain, bty="L")
#'

plot_data_assigned <- function(chain, legend=TRUE, col_early="hotpink",
                               col_late="skyblue3",
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
  xmax <- which(avg>0.999)[1]

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
