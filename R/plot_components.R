#' Plot number of data items assigned to each component for each posterior sample.
#'
#' @param chain A hdpSampleChain object including output from \code{\link{hdp_extract_components}}
#' @param legend Logical - should a legend be included? (default TRUE)
#' @param col_early Color ramp side for early posterior samples
#' @param col_late Color ramp side for late posterior samples
#' @param ... Other arguments to plot
#'  @seealso \code{\link{hdp_extract_components}},
#'  \code{\link{plot_comp_distn_bar}}, \code{\link{plot_dp_comp_exposure}}
#' @export
#' @examples
#' tcga_example_comp <- hdp_extract_components(tcga_example_chain)
#' plot_comp_size(tcga_example_comp, bty="L")

plot_comp_size <- function(chain, legend=TRUE, col_early="hotpink",
                           col_late="skyblue3", ...){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (length(comp_categ_counts(chain)) == 0) {
    stop("No component info for chain. First run hdp_extract_components")
  }
  if(class(legend) != "logical") stop("legend must be TRUE or FALSE")


  sums <- t(sapply(chain@comp_categ_counts, rowSums))

  # colour ramp across posterior samples
  cols <- colorRampPalette(colors=c(col_early, col_late))

  matplot(x=0:(nrow(sums)-1), sums, pch=1, col=cols(hdp_settings(chain)$n),
          xlab="Component", ylab="Number of data items", ...)

  if (legend) {
    legend("topright", col=c(col_early, col_late), pch=1,
           legend=c("Early samples", "Late samples"), bty="n")
  }
}



#' Barplot of the mean distribution over data categories for each component
#'
#' @param chain A hdpSampleChain object including output from \code{\link{hdp_extract_components}}
#' @param comp (Optional) Number(s) of the component(s) to plot (from 0 to the max component number).
#'  The default is to plot all components.
#' @param cat_names (Optional) Data category names to label the horizontal axis
#' @param grouping (Optional) A factor indicating data category groups.
#' @param col Either a single colour for all data categories, or a vector of
#'  colours for each group (in the same order as the levels of the grouping factor)
#' @param col_nonsig (Optional) Colour for any data category whose 95% credibility interval
#'  overlaps with zero (if set, overrides col argument)
#' @param show_group_labels Logical - should group labels be added to the top
#'  horizontal axis? (default FALSE) (only works if categories alreayd come in orders)
#' @param cred_int Logical - should 95\% credibility intervals be plotted? (default TRUE)
#' @param weights (Optional) Weights over the data categories to adjust their
#'  relative contribution (multiplicative)
#'  @seealso \code{\link{hdp_extract_components}}, \code{\link{plot_comp_size}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @export
#' @examples
#' tcga_example_comp <- hdp_extract_components(tcga_example_chain)
#' bases <- c('A', 'C', 'G', 'T')
#' trinuc_context <- paste0(rep(rep(bases, times=4), each=6),
#'                          rep(c('C', 'T'), each=48),
#'                          rep(bases, times=24))
#' group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
#'                            each=16))
#' plot_comp_distn_bar(tcga_example_comp, cat_names=trinuc_context,
#'                 grouping=group_factor, col=RColorBrewer::brewer.pal(6, "Set2"),
#'                col_nonsig="grey80", show_group_labels=TRUE)

plot_comp_distn_bar <- function(chain, comp=NULL, cat_names=NULL,
                            grouping=NULL, col="grey70", col_nonsig=NULL,
                            show_group_labels=FALSE, cred_int=TRUE,
                            weights=NULL){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (length(comp_categ_counts(chain)) == 0) {
    stop("No component info for chain. First run hdp_extract_components")
  }
  ncat <- numcateg(final_hdpState(chain))
  nsamp <- hdp_settings(chain)$n
  comp_distn <- comp_categ_distn(chain)
  ncomp <- nrow(comp_distn$mean)-1
  if (class(comp) != "NULL") {
    if (class(comp) != "numeric" | any(comp %% 1 != 0) |
          any(comp < 0) | any(comp > ncomp)) {
      stop(paste("comp must be integer between 0 and", ncomp))
    }
  }
  if (!class(cat_names) %in% c("character", "NULL") |
       !length(cat_names) %in% c(ncat, 0)) {
    stop("cat_names must be a character vector with one value for every
         data category, or NULL")
  }
  if (!class(grouping) %in% c("factor", "NULL") |
        !length(grouping) %in% c(ncat, 0)) {
    stop("grouping must be a factor with one value for every
         data category, or NULL")
  }
  if(class(show_group_labels) != "logical") {
    stop("show_group_labels must be TRUE or FALSE")
  }
  if(class(cred_int) != "logical") stop("cred_int must be TRUE or FALSE")
  if(!class(weights) %in% c("numeric", "NULL") |
       !length(weights) %in% c(ncat, 0)) {
    stop("weights must be a numeric vector with one value for every
         data category, or NULL")
  }

  # which components to plot
  if (is.null(comp)){
    comp_to_plot <- 0:ncomp
  } else {
    comp_to_plot <- comp
  }

  # colours for each category
  if (is.null(grouping)){
    cat_cols <- rep(col, ncat)
  } else {
    cat_cols <- col[grouping]
  }

  for (ii in comp_to_plot){
    # mean categorical distribution (sig), and credibility interval
    sig <- comp_distn$mean[as.character(ii),]
    ci <- comp_distn$cred.int[[as.character(ii)]]

    # adjust categories by weights if specified (lose credibility intervals at the moment)
    if(!is.null(weights)){
      sig <- sig %*% diag(weights)
      denom <- sum(sig)
      sig <- sig/denom
      ci <- NULL # not sure how to get cred int if adjusting with weights
    }

    sig <- as.vector(sig)

    # set categories whose credibility intervals hit zero to a different colour
    cat_cols_copy <- cat_cols
    if (!is.null(col_nonsig)){
      cat_cols_copy[which(ci[1,]==0)] <- col_nonsig
    }

    # max plotting height
    plottop <- max(0.2, ceiling(max(ci)/0.1)*0.1)

    # main barplot
    b <- barplot(sig, col=cat_cols_copy, xaxt='n', ylim=c(0,plottop*1.1),
                 border=NA, names.arg=rep('',96), xpd=F,
                 main=paste("Component", ii))

    # add credibility intervals
    if (cred_int & !is.null(ci)){
      segments(x0=b, y0=ci[1,], y1=ci[2,], col='grey30')
    }

    # add category names
    if (!is.null(cat_names)){
      mtext(cat_names, side=1, las=2, at=b, cex=0.7,
            family='mono', col=cat_cols)
    }

    # add group labels
    if(show_group_labels){
      gl <- rle(as.vector(grouping))

      glends <- cumsum(gl$lengths)
      glstarts <- c(1, glends[-length(glends)]+1)
      glcol <- col[as.factor(gl$values)]

      segments(x0=b[glstarts], x1=b[glends],
               y0=plottop, col=glcol, lwd=10)

      text(b[floor((glends+glstarts)/2)], y=plottop*1.05,
           labels=gl$values)
    }
  }
}




#' Plot the mean distribution over components for each specified DP
#'
#' @param chain A hdpSampleChain object including output from \code{\link{hdp_extract_components}}
#' @param dpindices Indices of DP nodes to plot
#' @param col Colours of each component, from 0 to the max number
#' @param dpnames (Optional) Names of the DP nodes
#' @param main_text (Optional) Text at top of plot
#' @param incl_numdata_plot Logical - should an upper barplot indicating the number of
#'  data items per DP be included? (Default TRUE)
#' @param incl_nonsig Logical - should components whose credibility intervals include 0
#'  be included (per DP)? (Default TRUE)
#'  @seealso \code{\link{hdp_extract_components}}, \code{\link{plot_comp_size}},
#'  \code{\link{plot_comp_distn_bar}}
#' @export
#' @examples
#' tcga_example_comp <- hdp_extract_components(tcga_example_chain)
#' plot_dp_comp_exposure(tcga_example_comp, 4:30,
#'                       RColorBrewer::brewer.pal(9, "Set3"))
#' plot_dp_comp_exposure(tcga_example_comp, 4:30,
#'                       RColorBrewer::brewer.pal(9, "Set3"),
#'                       incl_numdata_plot=FALSE, incl_nonsig=FALSE)

plot_dp_comp_exposure <- function(chain, dpindices, col, dpnames=NULL,
                                main_text=NULL, incl_numdata_plot=TRUE,
                                incl_nonsig=TRUE){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (length(comp_categ_counts(chain)) == 0) {
    stop("No component info for chain. First run hdp_extract_components")
  }
  dp_distn <- comp_dp_distn(chain)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  if (!is.numeric(dpindices) | any(dpindices %% 1 != 0) |
        any(dpindices < 1) | any(dpindices > ndp)) {
    stop(paste("dpindices must be integers between 1 and", ndp))
  }
  if (!length(col) >= ncomp) {
    stop(paste("col must specify at least", ncomp, "colours"))
  }
  if (!class(dpnames) %in% c("character", "NULL") |
        !length(dpnames) %in% c(length(dpindices), 0)) {
    stop("dpnames must be a character vector with
         same length as dpindices, or NULL")
  }
  if (!class(main_text) %in% c("character", "NULL") |
        !length(main_text) %in% c(1, 0)) {
    stop("main_text must be a character string, or NULL")
  }
  if(class(incl_numdata_plot) != "logical") {
    stop("incl_numdata_plot must be TRUE or FALSE")
  }
  if(class(incl_nonsig) != "logical") stop("incl_nonsig must be TRUE or FALSE")


  # save pre-existing par conditions, and reset on exit
  par_old <- par(no.readonly=TRUE)
  on.exit(par(par_old), add=TRUE)

  # Number of data items per DP
  dps <- dp(final_hdpState(chain))[dpindices]
  numdata <- sapply(dps, function(x) x@numdata)
  dp_order <- order(numdata, decreasing=TRUE)

  # mean exposures
  exposures <- t(dp_distn$mean[dpindices,])

  # adjust by only considering significantly non-zero exposures
  if (!incl_nonsig){
    cis <- dp_distn$cred.int[dpindices]
    nonsig <- lapply(cis, function(x) which(x[1,]==0))
    for (i in 1:length(nonsig)){
      exposures[nonsig[[i]],i] <- 0
    }
    #exposures <- exposures %*% diag(1/colSums(exposures))
  }

  # which components to include in this plot
  inc <- which(rowSums(exposures, na.rm=T)>0)

  num_leg_col <- max(1, floor(sqrt(length(inc)))-1)

  if (incl_numdata_plot){
    par(mfrow=c(2, 1), mar=c(1, 4, 2, 0.5), oma=c(1.5, 1.5, 1, 1), cex.axis=0.7)

    barplot(numdata[dp_order], main=main_text, col='gray', space=0, border=NA,
            names.arg=dpnames, ylab='Number of data items', las=2, cex.names=0.6,
            legend.text=names(inc), args.legend=list(fill=col[inc],
                                                          bty='n',
                                                          title='Component',
                                                          ncol=num_leg_col))

    barplot(exposures[inc, dp_order], space=0, col=col[inc], border=NA,
            names.arg=dpnames, ylab='Component exposure', las=2, cex.names=0.6)
  } else {

    par(cex.axis=0.7)
    # don't understand why legend.text needs rev() here and not in above case,
    # but seems to work?
    barplot(exposures[inc, dp_order], space=0, col=col[inc], border=NA, ylim=c(0,1),
            xlim=c(0,length(dpindices)*1.3), names.arg=dpnames,
            ylab='Component exposure', las=2, cex.names=0.6,
            legend.text=rev(names(inc)), args.legend=list(fill=col[inc],
                                                          bty='n',
                                                          title='Component',
                                                          ncol=num_leg_col))
  }

}
