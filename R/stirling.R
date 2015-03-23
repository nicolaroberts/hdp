# func to calculate the unsigned stirling numbers of the first kind.
stirling <- function(nn){
  
  if (!exists('maxnn', where=parent.frame(2))){
    assign('maxnn', 1, envir=parent.frame(2))
    assign('allss', list(1), envir=parent.frame(2))
  }
  
  
  maxnn.local <- get('maxnn', envir=parent.frame(2))
  
  if (nn > maxnn.local){ #only calculate this if needed
    allss.local <- get('allss', envir=parent.frame(2))
    
    allss.local[(length(allss.local)+1):nn] <- 0
    
    for (mm in (maxnn.local+1):nn){
      allss.local[[mm]] <- c(allss.local[[mm-1]]*(mm-1), 0) + c(0, allss.local[[mm-1]])
      mss <- max(allss.local[[mm]])
      allss.local[[mm]] <- allss.local[[mm]]/mss
    }
    
    assign('maxnn', nn, envir=parent.frame(2))
    assign('allss', allss.local, envir=parent.frame(2))
  }
  
  assign('nn', nn, envir=parent.frame(2))
  ss <- evalq(allss[[nn]], envir=parent.frame(2))
  return(ss)
}