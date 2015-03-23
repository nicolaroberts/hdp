# func to return weights from a random stick breaking process
randstick <- function(alpha,numclass){
  zz <- c(rbeta(numclass-1, 1, alpha), 1) #proportion of stick broken off
  beta <- zz*cumprod(c(1, 1-zz[1:numclass-1])) #amount of stick remaining
  return(beta)
}