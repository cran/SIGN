#' GeneMatching is a function to remove uncommon genes between a list of expression matrices
#' @param ExpList List of expression matrices
#' @return List of expression matrices restricted to the common genes between them

GeneMatching <- function(ExpList){
 
  CommonGenes <- as.character(rownames(ExpList[[1]]))
  for(ListIter in 1:length(ExpList)){
   CommonGenes <- intersect(CommonGenes, as.character(rownames(ExpList[[ListIter]])))
  }
  
  ExpList_Intersected <- list()
  for(ListIter in 1:length(ExpList)){
   ExpMat <- ExpList[[ListIter]]
   MatchInd <- which(as.character(rownames(ExpMat)) %in% intersect(CommonGenes, as.character(rownames(ExpMat))))
   ExpMat <- ExpMat[MatchInd,]
   ExpMat <- ExpMat[sort(as.character(rownames(ExpMat)), index.return = T)[[2]],]
   ExpList_Intersected[[ListIter]] <- ExpMat
  }
 
  names(ExpList_Intersected) <- names(ExpList)
  return(ExpList_Intersected)
 
}
