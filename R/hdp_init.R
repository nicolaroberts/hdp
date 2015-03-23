#' Initialises a HDP list object.
#' 
#' @param ppindex Integer (or vector of integers) specifying the indices of the parental process(es) for the initial DPs. The 'top' DP will have parent process '0'. 
#' @param cpindex Integer (or vector of integers) specifying the indices of the concentration parameters for the initial DPs. 
#' @param hh Numeric vector with one entry for each category in the dataset; specifies the parameters of the base Dirichlet distribution. To specify the uniform Dirichlet, set hh=rep(1/number_of_categories, number_of_categories).
#' @param alphaa Numeric vector of shape hyperparameters for the gamma priors for DP concentration parameters. 
#' @param alphab Numeric vector of inverse scale hyperparameters for the gamma priors for DP concentration parameters. 
#' @return A list with the initial HDP structure. 
#' @export
#' 

hdp_init <- function(ppindex,cpindex,hh,alphaa,alphab){
  
  check_ppindex(length(ppindex),ppindex)
  check_cpindex(max(cpindex),cpindex)
  
  HELDOUT <- 0
  
  hdp <- list()
  
  hdp$numdp         <- length(ppindex) 
  hdp$numconparam   <- length(alphaa) 
  hdp$base          <- list() 
  hdp$base$hh       <- hh
  hdp$base$classqq  <- newclass(hdp$base$hh) 
  hdp$base$numclass <- 0L 
  hdp$conparam      <- vector('list',hdp$numconparam) 
  hdp$dp            <- vector('list',hdp$numdp) 
  hdp$dpstate       <- as.integer(HELDOUT*rep(1,hdp$numdp)) 
  hdp$ppindex       <- as.integer(ppindex) 
  hdp$cpindex       <- as.integer(cpindex) 
  hdp$ttindex       <- as.integer(rep(0,hdp$numdp)) 
  
  tt <- rep(0,hdp$numconparam)
  for (jj in 1:hdp$numdp){
    cp              <- cpindex[jj]
    tt[cp]          <- tt[cp] + 1
    hdp$ttindex[jj] <- as.integer(tt[cp])
  }

  for (jj in 1:hdp$numdp){
    hdp$dp[[jj]]$datacc  <- vector('integer')
    hdp$dp[[jj]]$classnd <- 0L
    hdp$dp[[jj]]$classnt <- 0L
    hdp$dp[[jj]]$beta    <- 1
    hdp$dp[[jj]]$alpha   <- vector('numeric')
    hdp$dp[[jj]]$numdata <- 0L
    hdp$dp[[jj]]$datass  <- vector('integer')
  }

  for (cp in 1:hdp$numconparam){
    hdp$conparam[[cp]]$alphaa  <- alphaa[cp]
    hdp$conparam[[cp]]$alphab  <- alphab[cp]
    hdp$conparam[[cp]]$numdp   <- as.integer(sum(cpindex==cp))
    hdp$conparam[[cp]]$alpha   <- hdp$conparam[[cp]]$alphaa/hdp$conparam[[cp]]$alphab
    hdp$conparam[[cp]]$totalnd <- as.integer(rep(0,hdp$conparam[[cp]]$numdp))
    hdp$conparam[[cp]]$totalnt <- as.integer(rep(0,hdp$conparam[[cp]]$numdp))
  }
  
  return(hdp)
}