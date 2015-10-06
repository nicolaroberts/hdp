hdp_extract_comp_single <- function(chain, prop.ex=0.97, cos.merge=0.90){

  set.seed(hdp_seed(chain), kind="Mersenne-Twister", normal.kind="Inversion")

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (class(prop.ex) != "numeric" | prop.ex >=1 | prop.ex <= 0) {
    stop("prop.ex must be between 0 and 1")
  }
  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }

  # number of categories, DPs, and samples
  ncat <- numcateg(final_hdpState(chain))
  ndp <- numdp(final_hdpState(chain))
  nsamp <- hdp_settings(chain)$n

  # function to merge cols of a matrix according to new (numeric) column labels
  merge_cols <- function(mx, labels){
    colnames(mx) <- labels
    num_rows <- nrow(mx)
    # init matrix of 0s with new column names
    ans <- matrix(0, nrow=num_rows, ncol=length(unique(labels)),
                  dimnames=list(c(1:num_rows), sort(as.numeric(unique(labels))))
                  )

    for (label in colnames(ans)){
      ans[, label] <- rowSums(matrix(mx[, which(colnames(mx) == label)],
                                    nrow=num_rows))
    }
    return(ans)
  }


  # Step (1)
  # Make each ccc (clust_categ_counts) and
  # cdc (clust_dp_counts) and
  # cdw (clust_dp_weights) matrix have the
  # same number of columns
  maxclust <- max(numcluster(chain))
  clust_label <- 1:maxclust

  ccc_1 <- lapply(clust_categ_counts(chain), function(x){
    ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })

  cdc_1 <- lapply(clust_dp_counts(chain), function(x){
    ans <- cbind(x, matrix(0, nrow=ndp, ncol=(maxclust-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })

  cdw_1 <- lapply(clust_dp_weights(chain), function(x){
    ans <- cbind(x, matrix(0, nrow=ndp, ncol=(maxclust-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })



  # Step (2)
  # Match up raw clusters (matrix columns) across posterior samples (columns not
  # guaranteed to keep same component through all samples)
  # K-centroids clustering of all raw clusters with cannot-link constraints
  # within each posterior sample
  ccc_unlist <- t(do.call(cbind, ccc_1))
  groupfactor <- rep(1:nsamp, each=maxclust)

  ccc_clust <- flexclust::kcca(ccc_unlist, maxclust,
                           group=groupfactor,
                           family=flexclust::kccaFamily("kmeans",
                                             groupFun="differentClusters"))
  ccc_label <- split(flexclust::clusters(ccc_clust), groupfactor)

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

  cdw_2 <- mapply(function(cdw, label) {
                    colnames(cdw) <- label
                    cdw[, order(as.numeric(colnames(cdw)))]
                    },
                  cdw_1, ccc_label, SIMPLIFY=FALSE)

  remove(ccc_1, cdc_1, cdw_1, ccc_unlist, groupfactor, ccc_clust, ccc_label)

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
  }
  ccc_3 <- lapply(ccc_2, merge_cols, clust_label)

  remove(avgdistn, i, distns, clust_cos, clust_same, same)


  # Step (4)
  # Keep the largest components that explain at LEAST prop.ex
  # proportion of the data
  # and assign the leftover clusters to component '0'
  avg_prop_data <- rowMeans(sapply(ccc_3, colSums))/sum(ccc_3[[1]])
  cl_ordered <- names(avg_prop_data)[order(avg_prop_data, decreasing=T)]
  cl_below_prop <- which(cumsum(avg_prop_data[cl_ordered]) < prop.ex)
  if (length(cl_below_prop) == 0) {
    use_clust <- cl_ordered
  } else {
    use_clust <- cl_ordered[1:(max(cl_below_prop)+1)]
  }
  # update clust_label vector to reflect the merging of the small
  # clusters into the 0-th component
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_4 <- lapply(ccc_2, merge_cols, clust_label)
  cdc_4 <- lapply(cdc_2, merge_cols, clust_label)
  cdw_4 <- lapply(cdw_2, merge_cols, clust_label)

  remove(ccc_2, ccc_3, cdc_2, cdw_2, avg_prop_data,
         cl_below_prop, cl_ordered, use_clust, clust_label)


  # Step (5)
  # Rename overall component, order by number of data items (on average)
  # 0th component still goes first
  avg_ndi <- rowMeans(sapply(ccc_4, colSums))
  colorder <- c(1, setdiff(order(avg_ndi, decreasing=T), 1))

  ccc_5 <- lapply(ccc_4, function(x) {
    x <- x[, colorder]
    colnames(x) <- 0:(ncol(x)-1)
    return(x)
  })

  cdc_5 <- lapply(cdc_4, function(x) {
    x <- x[, colorder]
    colnames(x) <- 0:(ncol(x)-1)
    return(x)
  })

  cdw_5 <- lapply(cdw_4, function(x) {
    x <- x[, colorder]
    colnames(x) <- 0:(ncol(x)-1)
    return(x)
  })

  # number of components
  ncomp <- length(colorder)

  remove(ccc_4, cdc_4, cdw_4, avg_ndi, colorder)


  # Step (6)
  # Convert ccc into list of length ncomp, with matrices nsamp*ncat
  ccc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncat)), ncomp)
  for (i in 1:ncomp){
    ccc_ans[[i]] <- t(sapply(ccc_5, function(x) x[, i]))
  }

  # Convert cdc and cdw into list of length ndp, with matrices nsamp*ncomp
  cdc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncomp)), ndp)
  cdw_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncomp)), ndp)
  for (i in 1:ndp){
    cdc_ans[[i]] <- t(sapply(cdc_5, function(x) x[i, ]))
    cdw_ans[[i]] <- t(sapply(cdw_5, function(x) x[i, ]))
  }

  remove(ccc_5, cdc_5, cdw_5)


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

  # Step (8)
  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts and weights?)
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

  # add extracted components into chain hdpSampleChain slots
  chain@comp_settings <- list(prop.ex=prop.ex,
                              cos.merge=cos.merge)

  chain@comp_categ_counts <- ccc_ans
  chain@comp_dp_counts <- cdc_ans
  chain@comp_dp_weights <- cdw_ans
  chain@comp_categ_distn <- list(mean=ccc_mean,
                                 cred.int=ccc_credint)
  chain@comp_dp_distn <- list(mean=cdc_mean,
                              cred.int=cdc_credint)

  # check validity and return
  if (!validObject(chain)) warning("Not a valid hdpSampleChain object.")
  return(chain)
}
