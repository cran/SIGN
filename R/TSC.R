#' TSC is a function to calculate transcriprtional similarity coefficient between two biological pathways
#' @param PathwayExp1 Expression matrix of genes within the chosen pathway in the 1st set of samples
#' @param PathwayExp2 Expression matrix of genes within the chosen pathway in the 2nd set of samples
#' @return Transcriptional similarity coefficient
#' @examples
#' Pathway1_ExpMat <- matrix(runif(100,0,10), ncol = 10)
#' Pathway2_ExpMat <- matrix(runif(100,0,10), ncol = 10)
#' TSC(Pathway1_ExpMat, Pathway2_ExpMat)
#' @export
TSC <- function(PathwayExp1, PathwayExp2){
  AA <- PathwayExp1%*%t(PathwayExp1)
  BB <- PathwayExp2%*%t(PathwayExp2)
  AA0 <- AA - diag(diag(AA))
  BB0 <- BB - diag(diag(BB))
  TSC <- sum(diag(AA0%*%BB0))/sum(AA0^2)^.5/sum(BB0^2)^.5
  return(TSC)
}
