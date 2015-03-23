#' Posterior sampling across activated DPs. 
#' 
#' @param hdp A HDP list object. 
#' @param numburnin The number of burnin iterations.
#' @param numsample The number of posterior samples to collect.
#' @param numspace The number of iterations to leave between collected samples.
#' @param doconparam The number of iterations of concentration parameter sampling to perform between iterations (?)
#' @param dolik Leave as default (1)
#' @param dodebug Verbosity of debugging statements. 0 is least verbose, 4 is most verbose. 0 Highly recommended.
#' @return A list with seven elements: 'hdp' - the HDP list object at the end of the sampling process; 'lik' - a numeric vector of the HDP likelihood over the entire sampling process; 'numclass' - a numeric vector of the total number of classes/clusters in each posterior sample collected; 'classqq' - a list of matrices summarising the total number of data points of each category in each overall class/cluster at each posterior sample collected; 'classnd' - a list of matrices summarising the total number of data points in each DP in each overall class/cluster at each posterior sample collected; 'alpha'; and 'beta'.  
#' @export
hdp_posterior <- function(hdp,numburnin,numsample,numspace,doconparam,dolik=1,dodebug=0){
  
  mindifftime <- function(t1, t2){
    as.numeric(t2-t1, units='mins')
  }
  
  totiter <- numburnin + numsample*numspace
  lik     <- rep(0,numburnin+numsample*numspace)
  
  starttime <- Sys.time()

  # burn in
  output <- iterate(hdp,numburnin,doconparam,dolik,dodebug)
  hdp <- output[[1]]
  lik[1:numburnin] <- output[[2]]
  #report burn-in time
  prevtime <- Sys.time()
  print(sprintf("%d burn-in iterations in %1.2f mins", numburnin, mindifftime(starttime, prevtime)))
  curriter <- numburnin
  
  sample  <- rep(list(hdp_getstate(hdp)), numsample)
  
  for (samp in 1:numsample){
    output <- iterate(hdp,numspace,doconparam,dolik,dodebug)
    hdp <- output[[1]]
    lik[numburnin+(samp-1)*numspace+(1:numspace)] <- output[[2]]
    sample[[samp]] <- hdp_getstate(hdp)
    #report time every minute, in minute units. 
    tracktime <- Sys.time()
    curriter <- curriter + numspace
    if (mindifftime(prevtime,tracktime)>1){
      elapsedtime <- mindifftime(starttime, tracktime)
      print(sprintf("time %1.2f ETC %1.2f mins", elapsedtime, elapsedtime/curriter*totiter))
      prevtime <- tracktime
    }
  }
  
  numclass <- sapply(sample, function(x) x$numclass)
  classqq <- lapply(sample, function(x) x$classqq)
  classnd <- lapply(sample, function(x) t(sapply(x$classnd, unlist)))
  alpha <- sapply(sample, function(x) x$alpha)
  beta <- lapply(sample, function(x) x$beta)
  
  return(list(hdp=hdp, 
              lik=lik, 
              numclass=numclass, 
              classqq=classqq, 
              classnd=classnd, 
              alpha=alpha, 
              beta=beta))
}


