# func to initialise classqq matrix
# creates a 1-column matrix of zeros (one for every mutation type)
# hh is the vector of parameters for the base distribution
newclass <- function(hh){
  qq <- as.integer(rep(0, length(hh)))
  dim(qq) <- c(length(hh),1)
  return(qq)
}

# func to add data from one sample to the classqq matrix
# qq is the classqq matrix, counts the number of data items of
#   each mutation type in each class
# cc is a numeric vector of class allocations for each data point in one sample
# ss is a numeric vector of data observations (the mutation types) in one sample
additems <- function(qq,cc,ss){
  sc <- table(factor(ss, levels=1:nrow(qq)), factor(cc, levels=1:ncol(qq)))
  qq <- qq + sc
  qq <- matrix(qq, nrow=dim(qq)[1], ncol=dim(qq)[2])
  return(qq)
}

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
      for (jj in which(J == ii)){
        clik <- mm * weights[jj]
        clik <- cumsum(stirnum * exp(clik - max(clik)))
        numtable[jj] <- 1 + sum(runif(1) * clik[maxtable] > clik)
      }
    }
  }
  return(numtable)
}


# func to return weights from a random stick breaking process
randstick <- function(alpha,numclass){
  zz <- c(rbeta(numclass-1, 1, alpha), 1) #proportion of stick broken off
  beta <- zz * cumprod(c(1, 1 - zz[1:(numclass-1)])) #amount of stick remaining
  return(beta)
}

