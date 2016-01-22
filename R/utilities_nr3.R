# is any non-empty slot negative?
is.any.slot.negative <- function(object){
  slot_lengths <- sapply(slotNames(object), function(x) length(slot(object, x)))
  slots_not_empty <- which(slot_lengths >= 1)
  for (slotname in names(slots_not_empty)) {
    value <- slot(object, slotname)
    if (is.atomic(value)){
      if (any(value < 0, na.rm=TRUE)) return(TRUE)
    }
  }
  return(FALSE)
}


# function to merge cols of a matrix according to new (numeric) column labels
# ans will be sorted by numeric column name
merge_cols <- function(mx, labels){
  colnames(mx) <- labels
  num_rows <- nrow(mx)
  # init matrix of 0s with new column names
  ans <- matrix(0, nrow=num_rows, ncol=length(unique(labels)),
                dimnames=list(c(1:num_rows), sort(as.numeric(unique(labels))))
  )

  for (label in colnames(ans)){
    ans[, label] <- rowSums(matrix(mx[, which(colnames(mx) == label)],
                                   nrow=num_rows))
  }
  return(ans)
}


# function to merge elements of a list according to new (numeric) names
# ans will be sorted by numeric element names (is this still needed?)
merge_elems <- function(lst, labels){
  names(lst) <- labels
  # init list with new names
  ans <- vector("list", length(unique(labels)))
  names(ans) <- sort(as.numeric(unique(labels)))

  for (label in names(ans)){
    ans[[label]] <- Reduce("+", lst[which(labels == label)])
  }
  return(ans)
}
