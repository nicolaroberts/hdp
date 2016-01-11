# func to add new classes to a hdp data structure
qq_addclass <- function(hdp,newclass){
  HELDOUT <- 0L
  ACTIVE <- 2L

  oldcl <- hdp@base@numclass
  numcl <- oldcl + newclass
  hdp@base@numclass <- as.integer(numcl)
  hdp@base@classqq <- cbind(hdp@base@classqq, matrix(0L,
                                  nrow=nrow(hdp@base@classqq), ncol=newclass))

  for (jj in 1:hdp@numdp){
    if (hdp@dpstate[jj] != HELDOUT){
      hdp@dp[[jj]]@classnd[numcl+1] <- 0L
      hdp@dp[[jj]]@classnt[numcl+1] <- 0L
    }
    if (hdp@dpstate[jj] == ACTIVE){
      hdp@dp[[jj]]@beta[(oldcl+1):(numcl+1)] <- hdp@dp[[jj]]@beta[oldcl+1] *
        randstick(hdp@dp[[jj]]@alpha,newclass + 1)
    }
  }
  return(hdp)
}
