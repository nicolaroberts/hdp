# func to randomly assign a number of tables 
randnumtable <- function(weights,maxtable){
  numtable <- rep(0, length(maxtable))
  B <- unique(sort(maxtable))
  J <- match(maxtable, B)
  
  weights <- log(weights)
  
  for (ii in 1:length(B)){
    maxtable <- B[ii]
    if (maxtable > 0){
      mm <- 1:maxtable
      stirnum <- stirling(maxtable)
      for (jj in which(J==ii)){
        clik <- mm * weights[jj]
        clik <- cumsum(stirnum * exp(clik-max(clik)))
        numtable[jj] <- 1+sum(runif(1)*clik[maxtable] > clik)
      }
    }
  }
  return(numtable)
}