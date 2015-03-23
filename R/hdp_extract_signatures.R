#' Consolidates the large number of classes into a smaller number of overall signatures. 
#' 
#' @param output The output list object from hdp_posterior
#' @param prop.explained The proportion of data points to assign to overall signatures (e.g. 0.95)
#' @param cos.merge Merge signatures with cosine similarity greater than this (e.g. 0.90)
#' @return A list with three elements:  map_classes_to_sigs - the mapping from class columns in the original classqq to the signature columns in sigs_qq; sigs_qq - a list of matrices summarising the total number of data points of each category in each overall signature at each posterior sample collected (signature 0 includes all unassigned data points); and sigs_nd_by_dp - a list of matrices summarising the total number of data points assigned to each overall signature in each posterior sample collected, for each DP in the HDP (signature 0 includes all unassigned data points).  
#' @export

hdp_extract_signatures <- function(output, prop.explained, cos.merge){
  # number of categories and number of DPs
  num_cats <- nrow(output$classqq[[1]])
  num_dps <- nrow(output$classnd[[1]])
  num_samples <- length(output$classnd)
  
  #function to merge columns of a matrix according to column labels  
  merge_cols <- function(mx, labels){
    colnames(mx) <- labels
    num_rows <- nrow(mx)
    ans <- matrix(0, nrow=num_rows, ncol=length(unique(labels)), 
                  dimnames=list(c(1:num_rows), sort(as.numeric(unique(labels)))))
    
    for (label in colnames(ans)){
      ans[,label] <- rowSums(matrix(mx[,which(colnames(mx)==label)], nrow=num_rows))
    }
    return(ans)
  }
  
  # Step (1)
  # Make each classqq matrix have the same number of columns (end two cols are superflous)
  maxclass <- max(output$numclass)+1
  class_label <- 1:maxclass
  classqq_1 <- lapply(output$classqq, function(x) cbind(x, matrix(0, nrow=num_cats, ncol=(maxclass-ncol(x)))))
  
  # Step (2)
  # Merge the classqq columns with high cosine similarity. 
  avgdistn <- matrix(0, nrow=num_cats, ncol=maxclass)
  for (i in 1:maxclass){
    distns <- sapply(classqq_1, function(x) x[,i]/sum(x[,i]))
    avgdistn[,i] <- rowMeans(distns, na.rm=T)
  }
  class_cos <- lsa::cosine(avgdistn)
  qq_same <- (class_cos > cos.merge & lower.tri(class_cos))
  same <- which(qq_same, arr.ind=TRUE) # merge these columns
  # update class_label vector to reflect the merging of columns. 
  if (length(same)>0){
    for (i in 1:nrow(same)){
      class_label[same[i,1]] <- class_label[same[i,2]]
    }  
  }
  classqq_2 <- lapply(classqq_1, merge_cols, class_label)
  
  # Step (3)
  # Keep the largest signatures that explain upto prop.explained proportion of the data
  # and assign the leftover classes to signature '0'
  avg_prop_data <- rowMeans(sapply(classqq_2, colSums))/sum(classqq_2[[1]])
  classes_ordered <- names(avg_prop_data)[order(avg_prop_data, decreasing=T)]
  use_classes <- classes_ordered[which(cumsum(avg_prop_data[classes_ordered]) < prop.explained)]
  # update class_label vector to reflect the merging of the small classes into the 0-th signature
  class_label[which(!class_label %in% use_classes)] <- 0
  classqq_3 <- lapply(classqq_1, merge_cols, class_label)
  
  # Step (4)
  # Rename overall signatures
  rename <- function(matrix){
    matrix <- matrix[,order(as.numeric(colnames(matrix)))]
    colnames(matrix) <- 0:(ncol(matrix)-1)
    return(matrix)
  }
  classqq_ans <- lapply(classqq_3, rename)
  
  transformation <- sort(unique(class_label)) #old names
  names(transformation) <- as.numeric(colnames(classqq_ans[[1]])) #new names
  class_label_ans <- class_label
  for (j in 1:length(transformation)){
    class_label_ans[which(class_label==transformation[j])] <- names(transformation)[j]
  }
  
  # Step (5)
  # For each DP in each sample, get exposures to each global sig.
  num_sigs <- length(transformation) # number of sigs + 1 for the 0-th
  # Make each classnd matrix have the same number of columns
  classnd_1 <- lapply(output$classnd, function(x) cbind(x, matrix(0, nrow=num_dps, ncol=(maxclass-ncol(x)))))
  # Merge columns according to the new class_labels (signatures)
  classnd_2 <- lapply(classnd_1, merge_cols, class_label_ans)
  # Initilise list of matrices (one for each DP), each will contain mutation counts in each global sig for each posterior sample collected. 
  dp_classnd <- rep(list(matrix(0, nrow=num_samples, ncol=num_sigs)), num_dps)
  for (i in 1:num_dps){
    dp_classnd[[i]] <- t(sapply(classnd_2, function(x) x[i,]))
  }
  
  return(list(map_classes_to_sigs=class_label_ans, sigs_qq=classqq_ans, sigs_nd_by_dp=dp_classnd))  
}
