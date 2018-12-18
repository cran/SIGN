#' Genes_SimCal is a function to calculate similarity between a set of samples
#'  and 2 reference groups of samples
#' @param ExpMat_Test Expression matrix for the test samples for which SIGN will 
#' indetify the similarity with the 2 reference sataset
#' @param ExpMat_Ref1 Expression matrix for the 1st reference set fo samples
#' @param ExpMat_Ref2 Expression matrix for the 2nd reference set fo samples
#' @param RefIDs Annotations corresponding to the 2 expression matrices 
#' (1st and 2nd names are associated with the 1st and 2nd expression matrix and )
#' @param TestClassIter Index to be matched with RefIDs for removal of test samples
#'  from reference expression matrices
#' @param SampleIter Index of samples in the test expression matrix exist in 
#' referencece expression matrix 1 or 2
#' @return Vector of similarity between the target samples and the 2 reference sets

Genes_SimCal <- function(ExpMat_Test, ExpMat_Ref1, ExpMat_Ref2, RefIDs, TestClassIter, SampleIter){
 
  rownames(ExpMat_Ref1) <- NULL
  colnames(ExpMat_Ref1) <- NULL
  if(TestClassIter == RefIDs[1]){
   ExpMat_Ref1 <- ExpMat_Ref1[,-SampleIter]
  }
  
  rownames(ExpMat_Ref2) <- NULL
  colnames(ExpMat_Ref2) <- NULL
  if(TestClassIter == RefIDs[2]){
   ExpMat_Ref2 <- ExpMat_Ref2[,-SampleIter]
  }
  
  GeneVarSortInd <- sort(apply(cbind(ExpMat_Ref1, ExpMat_Ref2), 1, function(X){
    mad(na.omit(as.numeric(X)))}), decreasing = T, index.return = T)[[2]]
  
  BubbleSort_Vec1 <- c()
  BubbleSort_Vec2 <- c()

  Pearson_Vec1 <- c()
  Pearson_Vec2 <- c()

  NameVec <- c()
  for(GeneNum in c(1e3,nrow(ExpMat_Ref1))){
   NameVec <- c(NameVec, paste(GeneNum,"genes", sep = "", collapse = "")) 

   MostVarInd <- GeneVarSortInd[1:GeneNum]
   
   TargetGenes_RefMat1 <- ExpMat_Ref1[MostVarInd,]
   TargetGenes_RefMat2 <- ExpMat_Ref2[MostVarInd,]
   
   TargetGenes_TestMat <- ExpMat_Test[MostVarInd,]
   
   TargetGenes_Ref1Vec <- as.numeric(apply(TargetGenes_RefMat1,1,function(X){median(na.omit(as.numeric(X)))}))
   TargetGenes_Ref2Vec <- as.numeric(apply(TargetGenes_RefMat2,1,function(X){median(na.omit(as.numeric(X)))}))
   TargetGenes_TesTVec <- as.numeric(apply(TargetGenes_TestMat,1,function(X){median(na.omit(as.numeric(X)))}))

   BubbleSort_Vec1 <- c(BubbleSort_Vec1, as.numeric(BubbleSort(TargetGenes_Ref1Vec, TargetGenes_TesTVec)))
   Pearson_Vec1 <- c(Pearson_Vec1, as.numeric(cor.test(TargetGenes_Ref1Vec, TargetGenes_TesTVec)$estimate))
  
   BubbleSort_Vec2 <- c(BubbleSort_Vec2, as.numeric(BubbleSort(TargetGenes_Ref2Vec, TargetGenes_TesTVec)))
   Pearson_Vec2 <- c(Pearson_Vec2, as.numeric(cor.test(TargetGenes_Ref2Vec, TargetGenes_TesTVec)$estimate))
  }

  GeneSim_Out <- c(BubbleSort_Vec1, BubbleSort_Vec2, Pearson_Vec1, Pearson_Vec2)
  names(GeneSim_Out) <- c(paste(rep("BubbleSort1"), NameVec, sep = "_"),
                          paste(rep("BubbleSort2"), NameVec, sep = "_"),
                          paste(rep("Pearson1"), NameVec, sep = "_"),
                          paste(rep("Pearson2"), NameVec, sep = "_"))

  return(GeneSim_Out)
}
