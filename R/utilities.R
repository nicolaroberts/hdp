# is any non-empty slot negative?
is.any.slot.negative <- function(object){
  slot_lengths <- sapply(slotNames(object), function(x) length(slot(object, x)))
  slots_not_empty <- which(slot_lengths >= 1)
  for (slotname in names(slots_not_empty)) {
    value <- slot(object, slotname)
    if (is.atomic(value)){
      if (any(value < 0)) return(TRUE)
    }
  }
  return(FALSE)
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
