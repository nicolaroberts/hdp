#' Plot extracted components
#'
#' @param hdpsample A hdpSampleChain or hdpSampleMulti object including output
#'  from \code{\link{hdp_extract_components}}
#' @name plotcomp
#' @importFrom graphics abline axis barplot matplot mtext par plot points segments text
#' @importFrom grDevices colorRampPalette
#' @examples
#' mut_example_multi <- hdp_extract_components(mut_example_multi)
#' plot_comp_size(mut_example_multi, bty="L")
#' bases <- c("A", "C", "G", "T")
#' trinuc_context <- paste0(rep(rep(bases, times=6), each=4),
#'                          rep(c("C", "T"), each=48),
#'                          rep(bases, times=24))
#' group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
#'                            each=16))
#' plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
#'                 grouping=group_factor, col=RColorBrewer::brewer.pal(6, "Set2"),
#'                 col_nonsig="grey80", show_group_labels=TRUE)
#' plot_dp_comp_exposure(mut_example_multi, 5:30,
#'                       RColorBrewer::brewer.pal(12, "Set3"))
#' plot_dp_comp_exposure(mut_example_multi, 5:30,
#'                       RColorBrewer::brewer.pal(12, "Set3"),
#'                       incl_numdata_plot=FALSE, incl_nonsig=FALSE)
NULL
#> NULL


#' @param legend Logical - should a legend be included? (default TRUE)
#' @param col_a Color ramp side for early posterior samples (if hdpSampleChain) or first chain (if hdpSampleMulti)
#' @param col_b Color ramp side for late posterior samples (if hdpSampleChain) or last chain (if hdpSampleMulti)
#' @param xlab Horizontal axis label
#' @param ylab Vertical axis label
#' @param ... Other arguments to plot
#' @export
#' @rdname plotcomp
plot_comp_size <- function(hdpsample, legend=TRUE, col_a="hotpink",
                           col_b="skyblue3", xlab="Component",
                           ylab="Number of data items", ...){

  # input checks
  if (!class(hdpsample) %in% c("hdpSampleChain", "hdpSampleMulti")) {
    stop("hdpsample must have class hdpSampleChain or hdpSampleMulti")
  }
  # if (!validObject(hdpsample)) stop("hdpsample not valid") # too slow on big objects
  if (length(comp_categ_counts(hdpsample)) == 0) {
    stop("No component info for hdpsample. First run hdp_extract_components")
  }
  if(class(legend) != "logical") stop("legend must be TRUE or FALSE")


  sums <- t(sapply(comp_categ_counts(hdpsample), rowSums))

  # colour ramp across posterior samples
  cols <- colorRampPalette(colors=c(col_a, col_b))

  # colour vector

  if (class(hdpsample) == "hdpSampleChain") {
    mycols <- cols(ncol(sums))
    legtext <- c("Early samples", "Late samples")
  } else if (class(hdpsample) == "hdpSampleMulti") {
    postsamps <- sapply(chains(hdpsample), function(x) length(numcluster(x)))
    mycols <- rep(cols(length(postsamps)), postsamps)
    legtext <- c("First chain", "Last chain")
  }

  plot(x=jitter(rep(0:(nrow(sums)-1), ncol(sums))), xlab=xlab, xaxt="n",
       y=as.vector(sums), pch=1, col=mycols, ylab=ylab, ...)
  axis(1, at=0:(nrow(sums)-1), labels=rownames(sums))

  if (legend) {
    legend("topright", col=c(col_a, col_b), pch=1, legend=legtext, bty="n")
  }
}

#' @param comp (Optional) Number(s) of the component(s) to plot (from 0 to the max component number).
#'  The default is to plot all components.
#' @param cat_names (Optional) Data category names to label the horizontal axis
#' @param grouping (Optional) A factor indicating data category groups.
#' @param col Either a single colour for all data categories, or a vector of
#'  colours for each group (in the same order as the levels of the grouping factor)
#' @param col_nonsig (Optional) Colour for any data category whose 95\% credibility interval
#'  overlaps with zero (if set, overrides col argument)
#' @param show_group_labels Logical - should group labels be added to the top
#'  horizontal axis? (default FALSE) (only works if categories alreayd come in orders)
#' @param cred_int Logical - should 95\% credibility intervals be plotted? (default TRUE)
#' @param weights (Optional) Weights over the data categories to adjust their
#'  relative contribution (multiplicative)
#' @param plot_title (Optional) Character vector of custom plot titles (one for each component plotted)
#' @param group_label_height Multiplicative factor from top of plot for group label placement
#' @param cex.cat Expansion factor for the (optional) cat_names
#' @export
#' @rdname plotcomp
plot_comp_distn <- function(hdpsample, comp=NULL, cat_names=NULL,
                            grouping=NULL, col="grey70", col_nonsig=NULL,
                            show_group_labels=FALSE, cred_int=TRUE,
                            weights=NULL, plot_title=NULL,
                            group_label_height=1.05, cex.cat=0.7, ...){

  # input checks
  if (!class(hdpsample) %in% c("hdpSampleChain", "hdpSampleMulti")) {
    stop("hdpsample must have class hdpSampleChain or hdpSampleMulti")
  }
  # if (!validObject(hdpsample)) stop("hdpsample not valid") # too slow on big objects
  if (length(comp_categ_counts(hdpsample)) == 0) {
    stop("No component info for hdpsample. First run hdp_extract_components")
  }

  if (class(hdpsample) == "hdpSampleChain") {
    ncat <- numcateg(final_hdpState(hdpsample))
  } else if (class(hdpsample) == "hdpSampleMulti") {
    ncat <- numcateg(final_hdpState(chains(hdpsample)[[1]]))
  }

  comp_distn <- comp_categ_distn(hdpsample)
  ncomp <- nrow(comp_distn$mean)-1
  if (class(comp) != "NULL") {
    if ( !is.numeric(comp) | any(comp %% 1 != 0) |
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
    comp_to_plot <- rownames(comp_distn$mean)
  } else {
    comp_to_plot <- rownames(comp_distn$mean)[comp+1]
  }

  if(!class(plot_title) %in% c("character", "NULL") |
       !length(plot_title) %in% c(length(comp_to_plot), 0)){
    stop("plot_title must be a character vector with one value for every
         component being plotted, or NULL")
  }

  # colours for each category
  if (is.null(grouping)){
    cat_cols <- rep(col, ncat)
  } else {
    cat_cols <- col[grouping]
  }

  # main titles
  if (is.null(plot_title)){
    plot_title <- paste("Component", comp_to_plot)
  }
  names(plot_title) <- comp_to_plot

  for (ii in seq_along(comp_to_plot)){

    cname <- comp_to_plot[ii]

    # mean categorical distribution (sig), and credibility interval
    sig <- comp_distn$mean[cname,]
    ci <- comp_distn$cred.int[[cname]]

    # adjust categories by weights if specified (lose cred intervals though)
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
    plottop <- ceiling(max(ci)/0.1)*0.1

    # main barplot
    b <- barplot(sig, col=cat_cols_copy, xaxt="n", ylim=c(0,plottop*1.1),
                 border=NA, names.arg=rep("", ncat), xpd=F, las=1,
                 main=plot_title[ii], ...)

    # add credibility intervals
    if (cred_int & !is.null(ci)){
      segments(x0=b, y0=ci[1,], y1=ci[2,], col="grey30")
    }

    # add category names
    if (!is.null(cat_names)){
      mtext(cat_names, side=1, las=2, at=b, cex=cex.cat,
            family="mono", col=cat_cols)
    }

    # add group labels
    if(show_group_labels){
      gl <- rle(as.vector(grouping))

      glends <- cumsum(gl$lengths)
      glstarts <- c(1, glends[-length(glends)]+1)
      glcol <- col[as.factor(gl$values)]

      bd <- mean(diff(b))
      segments(x0=b[glstarts]-bd, x1=b[glends]+bd,
               y0=plottop, col=glcol, lwd=9, lend=1)

      text(b[floor(glends/2 + glstarts/2)], y=plottop*group_label_height,
           labels=gl$values)
    }
  }
}


#' @param dpindices Indices of DP nodes to plot
#' @param col_comp Colours of each component, from 0 to the max number
#' @param dpnames (Optional) Names of the DP nodes
#' @param main_text (Optional) Text at top of plot
#' @param incl_numdata_plot Logical - should an upper barplot indicating the number of
#'  data items per DP be included? (Default TRUE)
#' @param incl_nonsig Logical - should components whose credibility intervals include 0
#'  be included (per DP)? (Default TRUE)
#' @param incl_comp0 Logical - should component zero be plotted? (Default TRUE)
#' @param ylab_numdata Vertical axis label for numdata plot
#' @param ylab_exp Vertical exis label for exposure plot
#' @param leg.title Legend title
#' @param cex.names Expansion factor for bar labels (dpnames) in exposure plot
#' @param cex.axis Expansion factor for vertical-axis annotation
#' @param mar See ?par
#' @param oma See ?par
#' @export
#' @rdname plotcomp
plot_dp_comp_exposure <- function(hdpsample, dpindices, col_comp, dpnames=NULL,
                                main_text=NULL, incl_numdata_plot=TRUE,
                                incl_nonsig=TRUE, incl_comp0=TRUE,
                                ylab_numdata="Number of data items",
                                ylab_exp="Component exposure",
                                leg.title="Component", cex.names=0.6,
                                cex.axis=0.7, mar=c(1, 4, 2, 0.5),
                                oma=c(1.5, 1.5, 1, 1), ...){

  # input checks
  if (!class(hdpsample) %in% c("hdpSampleChain", "hdpSampleMulti")) {
    stop("hdpsample must have class hdpSampleChain or hdpSampleMulti")
  }
  # if (!validObject(hdpsample)) stop("hdpsample not valid") # too slow on big objects
  if (length(comp_categ_counts(hdpsample)) == 0) {
    stop("No component info for hdpsample. First run hdp_extract_components")
  }
  dp_distn <- comp_dp_distn(hdpsample)
  ndp <- nrow(dp_distn$mean)
  ncomp <- ncol(dp_distn$mean)
  if (!is.numeric(dpindices) | any(dpindices %% 1 != 0) |
        any(dpindices < 1) | any(dpindices > ndp)) {
    stop(paste("dpindices must be integers between 1 and", ndp))
  }
  if (!length(col_comp) >= ncomp) {
    stop(paste("col_comp must specify at least", ncomp, "colours"))
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
  if (class(hdpsample) == "hdpSampleChain") {
    dps <- dp(final_hdpState(hdpsample))[dpindices]
    pps <- ppindex(final_hdpState(hdpsample))[dpindices]
  } else if (class(hdpsample) == "hdpSampleMulti") {
    dps <- dp(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
    pps <- ppindex(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
  }

  numdata <- sapply(dps, function(x) x@numdata)
  dp_order <- order(numdata, decreasing=TRUE)

  # if incl_numdata_plot TRUE, throw error if one DP has no data associated
  if (incl_numdata_plot & any(numdata == 0)) {
    stop("Can't have incl_numdata_plot TRUE if
         one or more dpindices have no associated data item/s")
  }

  # if different parent indices, warning
  if (length(unique(pps)) > 1) {
    warning("some dpindices have different parent nodes,
            separate plots may be better")
  }

  # mean exposures
  exposures <- t(dp_distn$mean[dpindices,,drop=FALSE])
  # think I don't need this anymore with drop=FALSE above?
  # if (nrow(exposures)==1) {
  #   dim(exposures) <- rev(dim(exposures))
  #   rownames(exposures) <- colnames(dp_distn$mean)
  #   }

  # only include significantly non-zero exposures
  if (!incl_nonsig){
    cis <- dp_distn$cred.int[dpindices]
    nonsig <- lapply(cis, function(x) which(x[1,]==0))
    for (i in 1:length(nonsig)){
      exposures[nonsig[[i]],i] <- 0
    }
  }

  # exclude signature 0
  if (!incl_comp0){
    exposures <- exposures[-1,]
  }

  # which components to include in this plot
  inc <- which(rowSums(exposures, na.rm=T)>0)

  num_leg_col <- floor(sqrt(length(inc)))

  if (incl_numdata_plot){
    par(mfrow=c(2, 1), mar=mar, oma=oma, cex.axis=cex.axis, las=1)

    barplot(numdata[dp_order], main=main_text, col="gray", space=0, border=NA,
            names.arg="", ylab=ylab_numdata,
            legend.text=names(inc),
            args.legend=list(fill=col_comp[inc], bty="n", title=leg.title,
                             ncol=num_leg_col), ...)

    barplot(as.matrix(exposures[inc, dp_order, drop=FALSE]), space=0, col=col_comp[inc], border=NA,
            ylim=c(0, 1), names.arg=dpnames[dp_order], ylab=ylab_exp,
            cex.names=cex.names, ...)
  } else {

    par(cex.axis=cex.axis, las=2)
    # don't understand why legend.text needs rev() here and not in above case,
    # but seems to work?
    barplot(as.matrix(exposures[inc, dp_order, drop=FALSE]), space=0, col=col_comp[inc],
            border=NA, ylim=c(0, 1),
            xlim=c(0, length(dpindices) + num_leg_col + 1),
            names.arg=dpnames[dp_order],
            ylab=ylab_exp, cex.names=cex.names,
            legend.text=rev(names(inc)),
            args.legend=list(fill=col_comp[inc], bty="n", title=leg.title,
                             ncol=num_leg_col), ...)
  }

}
