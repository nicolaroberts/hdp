#' Extract major components from the raw clusters
#'
#' If prior components included via \code{\link{hdp_prior_init}} are
#' preserved by \code{hdp_extract_components}, they are prefixed with "P".
#' Any new components in this case are prefixed with "N".
#'
#' @param x hdpSampleChain or hdpSampleMulti object
#' @param cos.merge Merge components with cosine similarity above this threshold (default 0.90)
#' @param min.sample Components must have significant exposure in at least this many samples (i.e. those DP nodes with data assigned) (default 1)
#' @return A hdpSampleChain or hdpSampleMulti object updated with component information
#' @aliases hdp_extract_components
#' @seealso \code{\link{hdp_posterior}}, \code{\link{hdp_multi_chain}},
#'  \code{\link{plot_comp_size}}, \code{\link{plot_comp_distn}},
#'  \code{\link{plot_dp_comp_exposure}}
#' @import clue
#' @export
#' @examples
#' hdp_extract_components(mut_example_chain)
hdp_extract_components <- function(x, cos.merge=0.90, min.sample=1){

  # input checks
  if (class(x)=="hdpSampleChain") {
      warning('Extracting components on single posterior sampling chain. Recommend switching to multiple independent chains in a hdpSampleMulti object, see ?hdp_multi_chain')
      is_multi <- FALSE
    } else if (class(x)=="hdpSampleMulti") {
      is_multi <- TRUE
    } else {
      stop("x must have class hdpSampleChain or hdpSampleMulti")
    }

  if (!validObject(x)) stop("x not a valid object")

  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }
  if (min.sample %% 1 != 0 | min.sample < 1) {
    stop("min.sample must be a positive integer")
  }


  if (is_multi) {
    # list of hdpSampleChain objects
    chlist <- x@chains
    nch <- length(chlist)

    # set seed, get final state and number of posterior samples
    set.seed(sampling_seed(chlist[[1]]), kind="Mersenne-Twister", normal.kind="Inversion")
    finalstate <- final_hdpState(chlist[[1]])
    nsamp <- sum(sapply(chlist, function(x) hdp_settings(x)$n))

  } else {
    # set seed, get final state and number of posterior samples
    set.seed(sampling_seed(x), kind="Mersenne-Twister", normal.kind="Inversion")
    finalstate <- final_hdpState(x)
    nsamp <- hdp_settings(x)$n

  }

  # number of categories, DPs,data items at each DP, and frozen priors
  ncat <- numcateg(finalstate)
  ndp <- numdp(finalstate)
  numdata <- sapply(dp(finalstate), numdata)
  pseudo <- pseudoDP(finalstate)
  rm(finalstate)

  is_prior <- length(pseudo) > 0
  if (is_prior) {
    priorcc <- 1:length(pseudo)
  }


  # Step (1)
  # Make each ccc (clust_categ_counts) and
  # cdc (clust_dp_counts) matrix have the
  # same number of columns

  if(is_multi){

    maxclust <- max(sapply(chlist, function(x) max(numcluster(x))))
    clust_label <- 1:maxclust

    ccc_0 <- lapply(chlist, function(ch){
      lapply(clust_categ_counts(ch), function(x){
        ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
        return(ans[, -ncol(ans)])
      })
    })


    cdc_0 <- lapply(chlist, function(ch){
      lapply(clust_dp_counts(ch), function(x){
        ans <- cbind(x, matrix(0, nrow=ndp, ncol=(maxclust-ncol(x)+1)))
        return(ans[, -ncol(ans)])
      })
    })


    # if priors, remove pseudo-counts from ccc_0
    if (is_prior){
      pseudodata <- sapply(dp(final_hdpState(chlist[[1]]))[pseudo],
                           function(x) table(factor(x@datass, levels=1:ncat)))

      ccc_0 <- lapply(ccc_0, function(y) lapply(y, function(x) {
        x[,priorcc] <- x[,priorcc] - pseudodata
        return(x)
      }))
    }

    ccc_raw_avg_per_ch <- lapply(ccc_0, function(matlist){ Reduce('+', matlist)/length(matlist) })

    mclust <- ncol(ccc_raw_avg_per_ch[[1]])

    rapch_unlist <- t(do.call(cbind, ccc_raw_avg_per_ch))
    rapch_gf <- rep(1:nch, each=mclust)
    rapch_ic <- rep(1:mclust, times=nch)

    rapch_clust <- flexclust::kcca(rapch_unlist, k=rapch_ic,
                                   group=rapch_gf,
                                   family=flexclust::kccaFamily("kmedians",
                                                                groupFun="differentClusters"))

    rapch_label <- split(flexclust::clusters(rapch_clust), rapch_gf)

    ccc_1 <- Reduce('c', mapply(function(matlist, rank){
      lapply(matlist, function(mat){
        ans <- mat[,order(rank)]
        return(ans)
      })
    }, ccc_0, rapch_label, SIMPLIFY=FALSE))

    cdc_1 <- Reduce('c', mapply(function(matlist, rank){
      lapply(matlist, function(mat){
        ans <- mat[,order(rank)]
        return(ans)
      })
    }, cdc_0, rapch_label, SIMPLIFY=FALSE))

    remove(ccc_0, cdc_0, ccc_raw_avg_per_ch, rapch_unlist, rapch_gf, rapch_ic,
           rapch_clust, rapch_label, mclust)

  } else {

    maxclust <- max(numcluster(x))
    clust_label <- 1:maxclust

    ccc_1 <- lapply(clust_categ_counts(x), function(x){
      ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
      return(ans[, -ncol(ans)])
    })

    cdc_1 <- lapply(clust_dp_counts(x), function(x){
      ans <- cbind(x, matrix(0, nrow=ndp, ncol=(maxclust-ncol(x)+1)))
      return(ans[, -ncol(ans)])
    })

    # if priors, remove pseudo-counts from ccc_1
    if (is_prior){
      pseudodata <- sapply(dp(final_hdpState(x))[pseudo],
                           function(x) table(factor(x@datass, levels=1:ncat)))

      ccc_1 <- lapply(ccc_1, function(x) {
        x[,priorcc] <- x[,priorcc] - pseudodata
        return(x)
      })
    }

  }


  # Step (2)
  # Match up raw clusters (matrix columns) across posterior samples (columns not
  # guaranteed to keep same component through all samples)

  # K-centroids clustering of all raw clusters with cannot-link constraints
  # within each posterior sample, Manhattan distance and median centroid
  mclust <- ncol(ccc_1[[1]])

  if (mclust==1){
    ccc_label <- rep(1, length(ccc_1))

  } else{
    ccc_unlist <- t(do.call(cbind, ccc_1))
    groupfactor <- rep(1:nsamp, each=mclust)
    initial_clust <- rep(1:mclust, times=nsamp)

    ccc_clust <- flexclust::kcca(ccc_unlist, k=initial_clust,
                                 group=groupfactor,
                                 family=flexclust::kccaFamily("kmedians",
                                                              groupFun="differentClusters"))

    # want this plot to be as simple as possible
    # tmp <- matrix(flexclust::clusters(ccc_clust), byrow=T, ncol=mclust)
    # matplot(tmp, type='l', lty=1, main="kmedians")

    ccc_label <- split(flexclust::clusters(ccc_clust), groupfactor)

    remove(ccc_unlist, groupfactor, initial_clust, ccc_clust)

  }

  ccc_2 <- mapply(function(ccc, label) {
    colnames(ccc) <- label
    ccc[, order(as.numeric(colnames(ccc)))]
  },
  ccc_1, ccc_label, SIMPLIFY=FALSE)

  cdc_2 <- mapply(function(cdc, label) {
    colnames(cdc) <- label
    cdc[, order(as.numeric(colnames(cdc)))]
  },
  cdc_1, ccc_label, SIMPLIFY=FALSE)

  remove(ccc_1, cdc_1, ccc_label)

  # Step (3)
  # Merge the ccc columns with high cosine similarity.
  avgdistn <- matrix(0, nrow=ncat, ncol=maxclust)
  for (i in 1:maxclust){
    distns <- sapply(ccc_2, function(x) x[, i]/sum(x[, i]))
    avgdistn[, i] <- rowMeans(distns, na.rm=T)
  }
  clust_cos <- lsa::cosine(avgdistn)
  clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind=TRUE) # merge these columns

  # update clust_label vector to reflect the merging of columns.
  if (length(same)>0){
    for (i in 1:nrow(same)){
      clust_label[same[i, 1]] <- clust_label[same[i, 2]]
    }
    remove(i)
  }
  ccc_3 <- lapply(ccc_2, merge_cols, clust_label)
  cdc_3 <- lapply(cdc_2, merge_cols, clust_label)
  clust_label <- colnames(ccc_3[[1]])
  if (any(clust_label != colnames(cdc_3))) stop("problem in step 3!")

  remove(avgdistn, distns, clust_cos, clust_same, same, ccc_2, cdc_2)


  # Step (4)
  # Assign components with no *significantly* non-zero data categories
  # to component '0'
  use_clust <- c()
  for (ii in 1:ncol(ccc_3[[1]])) {
    compii <- sapply(ccc_3, function(x) x[,ii])
    lowerb <- apply(compii, 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, 0.95)[1], 3)
      }
    })
    if(any(lowerb>0)) use_clust <- c(use_clust, colnames(ccc_3[[1]])[ii])
  }

  # update clust_label vector
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_4 <- lapply(ccc_3, merge_cols, clust_label)
  cdc_4 <- lapply(cdc_3, merge_cols, clust_label)
  clust_label <- colnames(ccc_4[[1]])
  if (any(clust_label != colnames(cdc_4))) stop("problem in step 4!")

  remove(compii, ccc_3, cdc_3, ii, lowerb, use_clust)

  # Step (5)
  # Assign components with < min.sample *significantly* non-zero sample exposure
  # to component '0' (disregarding DP nodes with no data items (parent nodes))
  use_clust <- c()
  disregard <- which(numdata==0)
  for (ii in 1:ncol(cdc_4[[1]])) {
    compii <- sapply(cdc_4, function(x) x[,ii])
    lowerb <- apply(compii[-disregard,], 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, 0.95)[1], 3)
      }
    })
    if(sum(lowerb>0)>=min.sample) use_clust <- c(use_clust, colnames(cdc_4[[1]])[ii])
  }

  # update clust_label vector
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_5 <- lapply(ccc_4, merge_cols, clust_label)
  cdc_5 <- lapply(cdc_4, merge_cols, clust_label)
  clust_label <- colnames(ccc_5[[1]])
  if (any(clust_label != colnames(cdc_5))) stop("problem in step 5!")

  remove(compii, ccc_4, cdc_4, ii, lowerb, use_clust, disregard)

  # Step (6)
  # Rename overall component, order by number of data items (on average)
  # 0th component still goes first

  avg_ndi <- rowMeans(sapply(ccc_5, colSums))
  colorder <- c(1, setdiff(order(avg_ndi, decreasing=T), 1))


  # If priors,
  # update clust_label to reflect match (down to 0.85) with prior components
  if (is_prior) {

    nco <- length(clust_label)

    # average distribution over data categ for each component
    avgdistn <- sapply(2:nco, function(i){
      distns <- sapply(ccc_5, function(x) x[, i]/sum(x[, i]))
      ans <- rowMeans(distns, na.rm=T)
      return(ans)
    })

    # compare against original pseudodata distn
    match2_pseudo <- apply(avgdistn, 2, function(x) lsa::cosine(x, pseudodata))
    rownames(match2_pseudo) <- priorcc
    colnames(match2_pseudo) <- clust_label[-1]

    to_match <- TRUE
    while(to_match){
      if(!any(match2_pseudo>0.85)) {
        to_match <- FALSE
        break
      }

      best_match <- which(match2_pseudo==max(match2_pseudo), arr.ind = TRUE)
      old <- colnames(match2_pseudo)[best_match[2]]
      new <- paste0("P", rownames(match2_pseudo)[best_match[1]])
      clust_label[which(clust_label == old)] <- new

      match2_pseudo <- match2_pseudo[-best_match[1], -best_match[2], drop=FALSE]

    }
    suppressWarnings(rm(to_match, match2_pseudo, avgdistn, best_match, old, new, nco))

    ccc_5 <- lapply(ccc_5, function(x) {
      colnames(x) <- clust_label
      return(x)
    })

    cdc_5 <- lapply(cdc_5, function(x) {
      colnames(x) <- clust_label
      return(x)
    })

  }


  ccc_6 <- lapply(ccc_5, function(x) {
    x <- x[, colorder]
    if (is_prior){
      update <- setdiff(which(!grepl("P", colnames(x))), 1)
      if (length(update)>0){
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    } else {
      colnames(x) <- 0:(ncol(x)-1)
    }
    return(x)
  })

  cdc_6 <- lapply(cdc_5, function(x) {
    x <- x[, colorder]
    if (is_prior){
      update <- setdiff(which(!grepl("P", colnames(x))), 1)
      if (length(update)>0){
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    } else {
      colnames(x) <- 0:(ncol(x)-1)
    }
    return(x)
  })

  # number of components
  ncomp <- length(colorder)

  remove(ccc_5, cdc_5, avg_ndi, colorder)


  # Step (7)
  # Convert ccc into list of length ncomp, with matrices nsamp*ncat

  ccc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncat)), ncomp)
  for (i in 1:ncomp){
    ccc_ans[[i]] <- t(sapply(ccc_6, function(x) x[, i]))
  }
  names(ccc_ans) <- colnames(ccc_6[[1]])

  # Convert cdc into list of length ndp, with matrices nsamp*ncomp
  cdc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncomp)), ndp)
  for (i in 1:ndp){
    cdc_ans[[i]] <- t(sapply(cdc_6, function(x) x[i, ]))
  }

  remove(ccc_6, cdc_6)

  # Step (8)
  # Calculate mean and 95% credibility interval for each component's
  # categorical data distribution
  ccc_norm <- lapply(ccc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm=TRUE))

  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })

  # Step (8)
  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts)
  cdc_norm <- lapply(cdc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm=TRUE))

  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })

  # add extracted components into x hdpSampleChain slots
  x@numcomp <- as.integer(ncomp - 1)

  # proportion of data explained by extracted components?
  avcount <- colMeans(sapply(ccc_ans, rowSums, na.rm=TRUE), na.rm=TRUE)
  x@prop.ex <- round(1-avcount[1]/sum(avcount), 3)

  x@comp_cos_merge <- cos.merge

  x@comp_categ_counts <- ccc_ans
  x@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")
  x@comp_categ_distn <- list(mean=ccc_mean,
                                 cred.int=ccc_credint)
  x@comp_dp_distn <- list(mean=cdc_mean,
                              cred.int=cdc_credint)

  # check validity and return
  if (!validObject(x)) warning("Not a valid hdpSampleChain/Multi object.")
  return(x)
}
