#' ExpPhen_Matching is a function for matching samples between expression matrices 
#' and metadata matrix (clinical feature matrix)
#' @param ExpMat Matrix of expression of genes (samples in columns and genes in rows)
#' @param MetaMat Matrix of clinical features (samples in columns)
#' @param SamID_Meta Sample ID in MetaMat
#' @return List of expression matrix and metadata of the clinical information after 
#' matching patiend IDs between the expression and clinical information matrices

ExpPhen_Matching <- function(ExpMat, MetaMat, SamID_Meta){
 
 SamIntersect <- intersect(as.character(colnames(ExpMat)), as.character(MetaMat[,SamID_Meta]))

 ExpMat_Matched <- ExpMat[,which(as.character(colnames(ExpMat)) %in% SamIntersect)]
 MetaMat_Matched <- MetaMat[which(as.character(MetaMat[,SamID_Meta]) %in% SamIntersect),]

 ExpMat_Matched <- ExpMat_Matched[,sort(as.character(colnames(ExpMat_Matched)), index.return = T)[[2]]]
 MetaMat_Matched <- MetaMat_Matched[sort(as.character(MetaMat_Matched[,SamID_Meta]), index.return = T)[[2]],]

 ExpMeta_List <- list(ExpMat_Matched, MetaMat_Matched)
 names(ExpMeta_List) <- c("Expression", "Meta")

 return(ExpMeta_List)

}
