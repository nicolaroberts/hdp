# func to check ppindex (indicies of parental processes)
# must be non-negative integer, and a smaller index than the DP itself. 
check_ppindex <- function(numdp,ppindex){
  if (any(ppindex<0) | any(ppindex>=1:numdp) | any(ppindex!=ceiling(ppindex))){
    stop('Not valid ppindex')
  }
}

# func to check dpindex vector (indicies of dirichlet processes)
# must be non-negative integers no larger than the total number of DPs, with no duplicates
check_dpindex <- function(numdp,dpindex){
  if (any(dpindex<=0) | any(dpindex>numdp) | any(dpindex!=ceiling(dpindex)) | length(dpindex)!=length(unique(dpindex))){
    stop('Not valid dpindex')
  }
}

# func to check cpindex vector (indicies of concentration parameters)
# must be positive integers no larger than the total number of concentration parameters 
check_cpindex <- function(numconparam,cpindex){
  if (any(cpindex<1) | any(cpindex>numconparam) | any(cpindex!=ceiling(cpindex))) {
    stop('Not valid cpindex')
  }
}

# func to initialise classqq matrix
# creates a 1-column matrix of zeros (one for every mutation type)
# hh is the vector of parameters for the base distribution
newclass <- function(hh){
  qq <- as.integer(rep(0, length(hh)))
  dim(qq) <- c(length(hh),1)
  return(qq)
}

# func to add data from one sample to the classqq matrix
# qq is the classqq matrix, counts the number of data items of each mutation type in each class
# cc is a numeric vector of class allocations for each data point in one sample
# ss is a numeric vector of data observations (the mutation types) in one sample
additems <- function(qq,cc,ss){
  sc <- table(factor(ss, levels=1:nrow(qq)), factor(cc, levels=1:ncol(qq)))
  qq <- qq + sc
  return(qq)
}