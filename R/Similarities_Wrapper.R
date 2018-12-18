#' Similarities_Wrapper is wrapper to identify similarities between the expression of genes in target sample
#' and the reference expression matrix
#' @param ExpMat_Test Expression matrix of test samples
#' @param ExpMat_Ref Expression matrix of reference samples
#' @param GeneVec Vector of gene names
#' @param PathwaySet List of pathways containing gene annotations for each pathways
#' @param RefID Class of the reference set
#' @param TestClassIter Class of the test set 
#' (if it is the same as reference set, the target test sample will be removed fro the reference set)
#' @param SampleIter Target test sample in ExpMat_Testto be used for comparison with ExpMat_Ref
#' @return List of similarities between the target sample and the expression matrix of reference samples

Similarities_Wrapper <- function(ExpMat_Test, ExpMat_Ref, GeneVec, PathwaySet, RefID, TestClassIter, SampleIter){
 
  rownames(ExpMat_Ref) <- NULL

  if(TestClassIter == RefID){
   ExpMat_Ref <- ExpMat_Ref[,-SampleIter]
  }
  SimMat <- Pathway_similarity(matrix(ExpMat_Test[,SampleIter], ncol = 1), ExpMat_Ref, GeneVec, PathwaySet, RefID)

  GSVA_Out <- GSVA_Calculation(matrix(ExpMat_Test[,SampleIter], ncol = 1), ExpMat_Ref, GeneVec, PathwaySet, RefID)

  Similarities <- list(SimMat, GSVA_Out)
  names(Similarities) <- c("Pathway_SimMat", "GSVA_Out")
  return(Similarities)
}
