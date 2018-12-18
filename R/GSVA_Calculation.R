#' GSVA_Calculation is a function for Calculating correlation between expression 
#' level of pathways between 2 groups using GSVA
#' @param ExpMat1 Expression matrix of genes in the 1st group of sampls
#' @param ExpMat2 Expression matrix of genes in the 2nd group of sampls
#' @param GeneVec Name of genes in the same order as considered in ExpMat1 and ExpMat2
#' @param GeneSets List of genes within pathways
#' @param Name Name used for naming the columns of output matrix of correlation between the 2 groups
#' @return Similarity of the pathway between the two expression matrices based on 
#' pearson correlation, bubble sort, and wilcoxon paaired rank test using GSVA enrichment scores of pathways
#' @importFrom stats cor.test
#' @importFrom stats wilcox.test
#' @importFrom stats pchisq
#' @importFrom GSVA gsva

GSVA_Calculation <- function(ExpMat1, ExpMat2, GeneVec, GeneSets, Name="SampleComparison"){

 rownames(ExpMat1) <- GeneVec
 rownames(ExpMat2) <- GeneVec
 
 ExpMat_Combined <- cbind(ExpMat1, ExpMat2)
 GSVA_TMP <- invisible(gsva(ExpMat_Combined, GeneSets, method='gsva'))

 GSVAMat1 <- matrix(GSVA_TMP$es.obs[,1], ncol = 1)
 GSVAMat2 <- matrix(GSVA_TMP$es.obs[,2:ncol(ExpMat_Combined)], ncol = (ncol(ExpMat_Combined)-1)) 

 Pearson_Vec <- apply(GSVAMat1, 2, function(X){apply(GSVAMat2, 2, function(Y){as.numeric(cor.test(as.numeric(X),as.numeric(Y))$estimate)})})
 BubbleSort_Vec <- apply(GSVAMat1, 2, function(X){apply(GSVAMat2, 2, function(Y){as.numeric(BubbleSort(as.numeric(X),as.numeric(Y)))})})
 WilcoxonPaired_Vec <- apply(GSVAMat1, 2, function(X){apply(GSVAMat2, 2, function(Y){as.numeric(wilcox.test(as.numeric(X),as.numeric(Y), paired=TRUE)$p.value)})})

 GSVA_output <- c(median(Pearson_Vec), median(BubbleSort_Vec), (1 - pchisq( -2*sum(log(na.omit(WilcoxonPaired_Vec))), 2*length(na.omit(WilcoxonPaired_Vec)))))
 names(GSVA_output) <- paste(Name, c("GSVA_Pearson", "GSVA_BubbleSort", "GSVA_WilcoxonPaired"), sep = "_")

 return(GSVA_output)
}
