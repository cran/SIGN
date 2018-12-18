#' Pathway_similarity is a function for calculating correlation between expression level of pathways
#'  between 2 groups using all the available approaches in SIGN 
#' @param ExpMat1 Expression matrix of genes in the 1st group of sampls
#' @param ExpMat2 Expression matrix of genes in the 2nd group of sampls
#' @param GeneVec Name of genes in the same order as considered in ExpMat1 and ExpMat2
#' @param GeneSets List of genes within pathways
#' @param Name Name used for naming the columns of output matrix of correlation between the 2 groups
#' @return Similarity of the pathway between the two expression matrices using 
#' pearson correlation, bubble sort, and wilcoxon paaired rank test

Pathway_similarity <- function(ExpMat1, ExpMat2, GeneVec, GeneSets, Name){
  
  TSC_Vec <- c()
  BubbleSort_Vec <- c()
  Pearson_Vec <- c()
  for(PathwayIter in 1:length(GeneSets)){
   TargetGenes  <- GeneSets[[PathwayIter]]
   MatchedInd <- which(GeneVec %in% TargetGenes)

   if(length(MatchedInd) > 1){
    PathExpMat1 <- matrix(ExpMat1[MatchedInd,], ncol = ncol(ExpMat1))
    PathExpMat2 <- matrix(ExpMat2[MatchedInd,], ncol = ncol(ExpMat2))
    TSC_Vec <- c(TSC_Vec, TSC(PathExpMat1, PathExpMat2))

    MedianVecTMP1 <- as.numeric(apply(PathExpMat1,1 ,function(X){median(na.omit(as.numeric(X)))}))
    MedianVecTMP2 <- as.numeric(apply(PathExpMat2,1 ,function(X){median(na.omit(as.numeric(X)))}))

    BubbleSort_Vec <- c(BubbleSort_Vec, as.numeric(BubbleSort(MedianVecTMP1, MedianVecTMP2)))
    Pearson_Vec <- c(Pearson_Vec, as.numeric(cor.test(MedianVecTMP1, MedianVecTMP2)$estimate))
   }
 }

  Similarity_Output <- cbind(TSC_Vec, BubbleSort_Vec, Pearson_Vec)
  colnames(Similarity_Output) <- paste(Name, c("TSC_Vec", "BubbleSort_Vec", "Pearson_Vec"), sep = "_")

  return(Similarity_Output)
}
