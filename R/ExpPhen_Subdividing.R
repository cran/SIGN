#' ExpPhen_Subdividing is a function for grouping samples based on a clinical feature
#'  available in metadata matrix (clinical feature matrix)
#' @param ExpMeta_List List containing expression matrix and metadata matrix
#' @param SubDiv_ID Index of the target clinical feature in metadata matrix for samples grouping 
#' @return List of expression and clinical information of patients grouped based on the specified clinical feature

ExpPhen_Subdividing <- function(ExpMeta_List, SubDiv_ID){
  
  ExpMat <- ExpMeta_List[["Expression"]]
  MetaMat <- ExpMeta_List[["Meta"]]

  FeatureVec <- as.character(MetaMat[,SubDiv_ID])
  UniqueFeature <- unique(FeatureVec)
  
  ExpList_MatchFeat <- list()
  MetaList_MatchFeat <- list()
  for(FeatIter in 1:length(UniqueFeature)){
   MatchedInd <- which(FeatureVec == UniqueFeature[FeatIter])
   ExpList_MatchFeat[[FeatIter]] <- ExpMat[,MatchedInd]
   MetaList_MatchFeat[[FeatIter]] <- MetaMat[MatchedInd,]
  }
  names(ExpList_MatchFeat) <- UniqueFeature
  names(MetaList_MatchFeat) <- UniqueFeature

  ExpMeta_MatchFeat <- list(ExpList_MatchFeat, MetaList_MatchFeat)
  names(ExpMeta_MatchFeat) <- c(paste("ExpList_", SubDiv_ID, sep = ""),
                                paste("MetaList_", SubDiv_ID, sep = ""))

  return(ExpMeta_MatchFeat)
}
