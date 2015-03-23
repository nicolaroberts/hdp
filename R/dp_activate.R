#' Activates the DPs to be included in the posterior sampling process. 
#' @param hdp A HDP list object. 
#' @param dpindex A vector of indices indicating which DPs in the HDP object to activate.
#' @param initcc A positive integer specifying how many clusters to start with. 
#' @export
dp_activate <- function(hdp,dpindex,initcc){
  
  check_dpindex(hdp$numdp,dpindex)
  
  ACTIVE  <- 2L
  FROZEN  <- 1L
  HELDOUT <- 0L
  
  # initialize numclass and classqq
  if (initcc > hdp$base$numclass){
    hdp <- qq_addclass(hdp,initcc-hdp$base$numclass)
  }

  dpindex <- sort(dpindex)
  # initialize state of HDP
  for (kk in 1:length(dpindex)){
    jj <- dpindex[kk]
    pp <- hdp$ppindex[jj]
    cp <- hdp$cpindex[jj]
    tt <- hdp$ttindex[jj]
    
    if (hdp$dpstate[jj] == ACTIVE){
      stop('Trying to activate a DP that is already activated')
    }
    if (pp > 0){
      if (hdp$dpstate[pp] != ACTIVE){
        stop('Ancestors of to be activated DPs has to be already activated')
      }
    }

    if (hdp$dpstate[jj] == HELDOUT){
      hdp$dp[[jj]]$datacc <- as.integer(ceiling(runif(hdp$dp[[jj]]$numdata)*hdp$base$numclass))
      hdp$base$classqq <- additems(hdp$base$classqq, hdp$dp[[jj]]$datacc, hdp$dp[[jj]]$datass)
      hdp$dp[[jj]]$classnd <- tabulate(hdp$dp[[jj]]$datacc,nbins=hdp$base$numclass+1)
      hdp$dp[[jj]]$classnt <- 0L
    }
    
    if (pp == 0){
      hdp$dp[[jj]]$beta <- rep(1,hdp$base$numclass+1)/(hdp$base$numclass+1)
    } else {
      hdp$dp[[jj]]$beta <- hdp$dp[[pp]]$beta
    }
    hdp$dp[[jj]]$alpha <- hdp$conparam[[cp]]$alpha
    hdp$dpstate[jj]  <- ACTIVE
  }
    
  for (kk in length(dpindex):1){
    jj <- dpindex[kk]
    pp <- hdp$ppindex[jj]
    cp <- hdp$cpindex[jj]
    tt <- hdp$ttindex[jj]
    alpha <- hdp$dp[[jj]]$alpha
  
    if (pp == 0){
      hdp$dp[[jj]]$classnt <- as.integer(hdp$dp[[jj]]$classnd > 0)
    } else {
      hdp$dp[[pp]]$classnd <- as.integer(hdp$dp[[pp]]$classnd - hdp$dp[[jj]]$classnt)
      hdp$dp[[jj]]$classnt <- as.integer(randnumtable(alpha*hdp$dp[[pp]]$beta, hdp$dp[[jj]]$classnd))
      hdp$dp[[pp]]$classnd <- as.integer(hdp$dp[[pp]]$classnd + hdp$dp[[jj]]$classnt)
    }

    hdp$conparam[[cp]]$totalnd[tt] <- as.integer(sum(hdp$dp[[jj]]$classnd))
    hdp$conparam[[cp]]$totalnt[tt] <- as.integer(sum(hdp$dp[[jj]]$classnt))
  }
  return(hdp)
}

