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
                                    nrow=num_rows), na.rm=TRUE)
    }
    return(ans)
  }


  # function to merge elements of a list according to new (numeric) names
  merge_elems <- function(lst, labels){
    names(lst) <- labels
    # init list with new names
    ans <- vector("list", length(unique(labels)))
    names(ans) <- unique(labels)

    for (label in names(ans)){
      ans[[label]] <- Reduce("+", lst[which(labels == label)])
    }
    return(ans[order(as.numeric(names(ans)))])
  }

  chains@comp_settings <- list(cos.merge=cos.merge)

  # list of hdpSampleChain objects
  chlist <- chains@chains
  nch <- length(chlist)

  # extract components on each chain
  for (i in 1:nch){

    if (!redo & length(comp_categ_counts(chlist[[i]])) > 0) next

    chlist[[i]] <- hdp_extract_comp_single(chlist[[i]],
                                           prop.ex=0.999, cos.merge=cos.merge)
  }

  chains@chains <- chlist


  # Consolidate components across chains.
  comp_mean_dist <- lapply(chlist, function(x) t(comp_categ_distn(x)$mean))
  nraw_comp <- sapply(comp_mean_dist, ncol) - 1

  # compare against chain with minimum number of components
  # and match to closest if cos similarity > 0.9
  baseline <- which.min(nraw_comp)

  for (i in setdiff(1:nch, baseline)) {
    new_col_names <- rep("0", nraw_comp[i])
    for (j in 1:nraw_comp[i]) {
      cos_sim <- apply(comp_mean_dist[[baseline]], 2, lsa::cosine,
                       comp_mean_dist[[i]][, as.character(j)])
      if (max(cos_sim, na.rm=TRUE) > 0.9) {
        new_col_names[j] <- names(which.max(cos_sim))
      }
    }
    new_col_names <- c("0", new_col_names)
    colnames(comp_mean_dist[[i]]) <-  new_col_names
  }

  comp_mapping <- lapply(comp_mean_dist, colnames)

  remove(comp_mean_dist, new_col_names, baseline, cos_sim, i, j)

  # if any colname not present in all chains, set to zero.
  comp_count <- table(unlist(sapply(comp_mapping, unique)))
  discard <- which(comp_count < nch)

  comp_mapping <- lapply(comp_mapping, function(x) {
    replace(x, which(x %in% names(discard)), "0")
    })

  remove(comp_count, discard)



  # Now have comp_mapping, consolidate all stats from all chains

  ncomp <- length(unique(comp_mapping[[1]]))

  # comp_categ_counts (rename overall comps so sorted by mean size)
  ccclist <- lapply(chlist, comp_categ_counts)
  cccmerge <- mapply(merge_elems, ccclist, comp_mapping)
  ranks <- rank(rowMeans(matrix(sapply(cccmerge, sum), nrow=ncomp)))
  names(ranks) <- rownames(cccmerge)

  dummy <- 1
  while (any(diff(ranks)[-1] > 0)){
    dummy <- dummy + 1
    which.switch <- which(diff(ranks)[-1] > 0)[1]
    to.switch <- names(ranks)[which.switch+1:2]
    one <- lapply(comp_mapping, function(x) which(x == to.switch[1]))
    two <- lapply(comp_mapping, function(x) which(x == to.switch[2]))
    comp_mapping <- mapply(function(x, y) {
      x[y] <- to.switch[2]
      return(x)
      }, comp_mapping, one, SIMPLIFY=FALSE)

    comp_mapping <- mapply(function(x, y) {
      x[y] <- to.switch[1]
      return(x)
    }, comp_mapping, two, SIMPLIFY=FALSE)

    cccmerge <- mapply(merge_elems, ccclist, comp_mapping)
    ranks <- rank(rowMeans(matrix(sapply(cccmerge, sum), nrow=ncomp)))
    names(ranks) <- rownames(cccmerge)

    if (dummy > 10) break

  }

  ccc <- vector("list", nrow(cccmerge))
  names(ccc) <- rownames(cccmerge)
  for (i in 1:length(ccc)){
    ccc[[i]] <- Reduce(rbind, cccmerge[i, ])
  }
  names(ccc) <- 0:(ncomp-1)
  remove(i, ccclist, cccmerge)


  # comp_dp_counts
  cdclist <- lapply(chlist, comp_dp_counts)
  cdcmerge <- mapply(function(x, y) lapply(x, merge_cols, y), cdclist,
                     comp_mapping, SIMPLIFY=FALSE)
  cdc <- vector("list", length(cdcmerge[[1]]))
  for (i in 1:length(cdc)) {
    cdc[[i]] <- Reduce(rbind, lapply(cdcmerge, `[[`, i))
    colnames(cdc[[i]]) <- 0:(ncomp-1)
  }
  remove(i, cdclist, cdcmerge)

  # comp_dp_weights
  cdwlist <- lapply(chlist, comp_dp_weights)
  cdwmerge <- mapply(function(x, y) lapply(x, merge_cols, y), cdwlist,
                     comp_mapping, SIMPLIFY=FALSE)
  cdw <- vector("list", length(cdwmerge[[1]]))
  for (i in 1:length(cdw)) {
    cdw[[i]] <- Reduce(rbind, lapply(cdwmerge, `[[`, i))
    colnames(cdw[[i]]) <- 0:(ncomp-1)
  }
  remove(i, cdwlist, cdwmerge)

  # Calculate mean and 95% credibility interval for each component's
  # categorical data distribution
  ccc_norm <- lapply(ccc, function(x) x/rowSums(x, na.rm=TRUE))

  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm=TRUE))
  rownames(ccc_mean) <- 0:(ncomp-1)

  ccc_credint <- lapply(ccc_norm, function(x) {apply(x, 2, function(y) {
    samp <- coda::as.mcmc(y)
    if (sum(!is.nan(samp)) ==  1) {
      c(NaN, NaN)
    } else {
      round(coda::HPDinterval(samp), 3)
    }})})
  names(ccc_credint) <- 0:(ncomp-1)

  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts and weights?)
  cdc_norm <- lapply(cdc, function(x) x/rowSums(x, na.rm=TRUE))

  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm=TRUE))

  cdc_credint <- lapply(cdc_norm, function(x) {apply(x, 2, function(y) {
    samp <- coda::as.mcmc(y)
    if (sum(!is.nan(samp)) ==  1) {
      c(NaN, NaN)
    } else {
      round(coda::HPDinterval(samp), 3)
    }})})


  # proportion of data explained by extracted components?
  avcount <- colMeans(sapply(ccc, rowSums, na.rm=TRUE), na.rm=TRUE)
  chains@prop.ex <- round(1-avcount[1]/sum(avcount), 3)

  # add extracted components into hdpSampleMulti slots
  chains@comp_categ_counts <- ccc
  chains@comp_dp_counts <- cdc
  chains@comp_dp_weights <- cdw
  chains@comp_categ_distn <- list(mean=ccc_mean,
                                 cred.int=ccc_credint)
  chains@comp_dp_distn <- list(mean=cdc_mean,
                              cred.int=cdc_credint)

# check validity and return
  if (!validObject(chains)) warning("Not a valid hdpSampleMulti object.")
  return(chains)
}
