hdp_extract_comp_multi <- function(chains, cos.merge=0.90, redo=TRUE){

  # input checks
  if (class(chains) != "hdpSampleMulti") {
    stop("chains must have class hdpSampleMulti")
  }
  if (!validObject(chains)) stop("chains not a valid hdpSampleMulti object")
  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }
  if (class(redo) != "logical") stop("redo must be TRUE or FALSE")

  # list of hdpSampleChain objects
  chlist <- chains@chains
  nch <- length(chlist)

  # extract components within each chain
  for (i in 1:nch){
    if (!redo & length(comp_categ_counts(chlist[[i]])) > 0) next
    chlist[[i]] <- hdp_extract_comp_single(chlist[[i]], cos.merge=cos.merge)
  }
  chains@chains <- chlist

  # number of categories, DPs, samples, and data items at each DP
  ncat <- numcateg(final_hdpState(chlist[[1]]))
  ndp <- numdp(final_hdpState(chlist[[1]]))
  nsamp <- sapply(chlist, function(x) hdp_settings(x)$n)
  numdata <- sapply(dp(final_hdpState(chlist[[1]])), numdata)

  remove(i)

  # Step (1) Match up components across chains
  # Get average comp_categ_counts matrix for each chain (average over posterior samples)
  ccc_0 <- lapply(chlist, function(x) apply(simplify2array(comp_categ_counts(x)), 2:3, mean))

  # Make each ccc matrix have the same number of columns
  maxcomp <- max(sapply(ccc_0, ncol))
  ccc_1 <- lapply(ccc_0, function(x){
    ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxcomp-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })

  remove(ccc_0)

  # K-centroids clustering of all components across chains
  # with cannot-link constraints within each chain.
  # Manhattan distance and median centroid.
  # exclude component 0 from each chain before matching
  ccc_1_no0 <- lapply(ccc_1, function(x) x[,-1])
  ccc_unlist <- t(do.call(cbind, ccc_1_no0))
  groupfactor <- rep(1:nch, each=maxcomp-1)
  initial_clust <- rep(1:(maxcomp-1), times=nch)

  ccc_clust <- flexclust::kcca(ccc_unlist, k=initial_clust,
                               group=groupfactor,
                               family=flexclust::kccaFamily("kmedians",
                                                    groupFun="differentClusters"))

  # want this plot to be as simple as possible
  # tmp <- matrix(flexclust::clusters(ccc_clust), byrow=T, ncol=maxcomp-1)
  # matplot(tmp, type='l', lty=1, main="kmedians")

  comp_mapping <- split(flexclust::clusters(ccc_clust), groupfactor)
  comp_mapping <- lapply(comp_mapping, function(x) c(0, x))

  remove(ccc_1, ccc_1_no0, ccc_unlist, groupfactor, initial_clust, ccc_clust)

  # Step (2) Consolidate all ccc and cdc stats from all chains
  # Get comp_categ_counts list for each chain,
  # make the same length, and reorder elements so they match up to comp_mapping
  ccclist <- lapply(chlist, comp_categ_counts)

  cccmerge <- mapply(function(x, y, z){
    ans <- c(x, rep(list(matrix(0, nrow=y, ncol=ncat)), (maxcomp-length(x)+1)))
    ans <- ans[-length(ans)]
    names(ans) <- z
    ans <- ans[order(as.numeric(names(ans)))]
    return(ans)
  }, ccclist, nsamp, comp_mapping, SIMPLIFY = TRUE)

  # rbind matrices from different chains
  ccc_2 <- vector("list", nrow(cccmerge))
  names(ccc_2) <- rownames(cccmerge)
  for (i in 1:length(ccc_2)){
    ccc_2[[i]] <- Reduce(rbind, cccmerge[i, ])
  }

  # Get comp_dp_counts list for each chain,
  # make number of columns the same, and reorder columns to match comp_mapping
  cdclist <- lapply(chlist, comp_dp_counts)
  cdcmerge <- mapply(function(x, y, z) {
    x <- lapply(x, function(w){
            ans <- cbind(w, matrix(0, nrow=y, ncol=(maxcomp-ncol(w)+1)))
            ans <- ans[, -ncol(ans)]
            colnames(ans) <- z
            ans <- ans[,order(as.numeric(colnames(ans)))]
            return(ans)
          })
  }, cdclist, nsamp, comp_mapping, SIMPLIFY=FALSE)

  # rbind matrices from different chains
  cdc_2 <- vector("list", length(cdcmerge[[1]]))
  for (i in 1:length(cdc_2)) {
    cdc_2[[i]] <- Reduce(rbind, lapply(cdcmerge, `[[`, i))
    colnames(cdc_2[[i]]) <- 0:(maxcomp-1)
  }

  remove(i, cdclist, cdcmerge, cccmerge, ccclist, chlist, comp_mapping, maxcomp)

  comp_label <- names(ccc_2)
  if(any(comp_label != colnames(cdc_2[[1]]))) stop("problem in step 2!")

  # Step (3)
  # Merge components with high cosine similarity in ccc
  avgdistn <- sapply(ccc_2, function(x) colMeans(x)/sum(colMeans(x)))
  clust_cos <- lsa::cosine(avgdistn)
  clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind=TRUE) # merge these

  # how best to track and update components??
  # update clust_label vector to reflect the merging of columns.
  if (length(same)>0){
    for (i in 1:nrow(same)){
      comp_label[same[i, 1]] <- comp_label[same[i, 2]]
    }
    remove(i)
  }
  ccc_3 <- merge_elems(ccc_2, comp_label)
  cdc_3 <- lapply(cdc_2, merge_cols, comp_label)

  remove(ccc_2, cdc_2, avgdistn, clust_cos, clust_same, same)

  comp_label <- names(ccc_3)
  if (any(comp_label != colnames(cdc_3))) stop("problem in step 3!")


  # Step (4)
  # Assign components with no *significantly* non-zero data categories
  # to component '0'
  use_comp <- c()
  for (ii in 1:length(ccc_3)) {
    lowerb <- apply(ccc_3[[ii]], 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (sum(!is.nan(samp)) ==  1) {
        NaN
      } else {
        round(coda::HPDinterval(samp)[1], 3)
      }
    })
    if(any(lowerb>0)) use_comp <- c(use_comp, names(ccc_3)[ii])
  }
  # update clust_label vector
  comp_label[which(!comp_label %in% use_comp)] <- 0
  ccc_4 <- merge_elems(ccc_3, comp_label)
  cdc_4 <- lapply(cdc_3, merge_cols, comp_label)

  remove(ccc_3, cdc_3, ii, lowerb, use_comp)

  comp_label <- names(ccc_4)
  if (any(comp_label != colnames(cdc_4))) stop("problem in step 4!")

  # Step (5)
  # Assign components with no *significantly* non-zero sample exposure
  # to component '0' (disregarding DP nodes with no data items (parent nodes))
  use_comp <- c()
  disregard <- which(numdata==0)
  for (ii in 1:ncol(cdc_4[[1]])) {
    compii <- t(sapply(cdc_4, function(x) x[,ii]))
    lowerb <- apply(compii[-disregard,], 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (sum(!is.nan(samp)) ==  1) {
        NaN
      } else {
        round(coda::HPDinterval(samp)[1], 3)
      }
    })
    if(any(lowerb>0)) use_comp <- c(use_comp, colnames(cdc_4[[1]])[ii])
  }
  # update comp_label vector
  comp_label[which(!comp_label %in% use_comp)] <- 0
  ccc_5 <- merge_elems(ccc_4, comp_label)
  cdc_5 <- lapply(cdc_4, merge_cols, comp_label)

  remove(compii, ccc_4, cdc_4, ii, lowerb, use_comp)

  comp_label <- names(ccc_5)
  if (any(comp_label != colnames(cdc_5))) stop("problem in step 5!")

  # Step (6)
  # Rename overall component, order by number of data items (on average)
  # 0th component still goes first
  avg_ndi <- sapply(ccc_5, function(x) mean(rowSums(x)))
  comporder <- c(1, setdiff(order(avg_ndi, decreasing=T), 1))

  ccc_ans <- ccc_5[comporder]
  names(ccc_ans) <- 0:(length(ccc_ans)-1)

  cdc_ans <- lapply(cdc_5, function(x) {
    x <- x[, comporder]
    colnames(x) <- 0:(ncol(x)-1)
    return(x)
  })

  # number of components
  ncomp <- length(comporder)

  remove(ccc_5, cdc_5, avg_ndi, comporder)


  # Step (7)
  # Calculate mean and 95% credibility interval for each component's
  # categorical data distribution
  ccc_norm <- lapply(ccc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm=TRUE))
  rownames(ccc_mean) <- 0:(ncomp-1)

  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (sum(!is.nan(samp)) ==  1) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp), 3)
      }
    })
  })
  names(ccc_credint) <- 0:(ncomp-1)

  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts)
  cdc_norm <- lapply(cdc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm=TRUE))

  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (sum(!is.nan(samp)) ==  1) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp), 3)
      }
    })
  })

  # add extracted components into hdpSampleMulti slots
  chains@numcomp <- as.integer(ncomp - 1)

  # proportion of data explained by extracted components?
  avcount <- colMeans(sapply(ccc_ans, rowSums, na.rm=TRUE), na.rm=TRUE)
  chains@prop.ex <- round(1-avcount[1]/sum(avcount), 3)

  chains@comp_cos_merge <- cos.merge

  chains@comp_categ_counts <- ccc_ans
  chains@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")
  chains@comp_categ_distn <- list(mean=ccc_mean,
                                 cred.int=ccc_credint)
  chains@comp_dp_distn <- list(mean=cdc_mean,
                              cred.int=cdc_credint)

# check validity and return
  if (!validObject(chains)) warning("Not a valid hdpSampleMulti object.")
  return(chains)
}
