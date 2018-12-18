#' BubbleSort is a function for calculating bubble sort correlation between two vectors
#' @param Vec1 Vector of values of 1st feature across samples
#' @param Vec2 Vector of values of 2nd feature across samples
#' @return Bubble sort similarity between the two vectors
#' @importFrom stats cov
#' @example
#' Pathway1_ExpVec <- runif(10,0,10)
#' Pathway2_ExpVec <- runif(10,0,10)
#' BubbleSort(Pathway1_ExpVec, Pathway2_ExpVec)
BubbleSort <- function(Vec1,Vec2){
  if(length(Vec1)!=length(Vec2)){
    stop("permutations need to be of same length!")
  }
  
  Vec1_rank <- rank(Vec1)
  Vec2_rank <- rank(Vec2)

  return(((choose(length(Vec1_rank),2) - cov(Vec1_rank,Vec2_rank,method='kendall')/2)/2)/choose(length(Vec1_rank),2))
}
