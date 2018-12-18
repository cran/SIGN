#' EventRenaming is a function for changing annotation of censored samples to 0
#'  and dead samples to 1 for survival analysis
#' @param EventVec Status vector for all of the samples (patients) including both samples undergone an event or censored
#' @param Censored_Annot Index of samples censored in the dataset
#' @return Vector of events including 0 for censoring and 1 for death

EventRenaming <- function(EventVec, Censored_Annot){

 CensoredInd <- which(EventVec == Censored_Annot)

 EventVec_Renamed <- EventVec
 EventVec_Renamed[CensoredInd] <- 0
 EventVec_Renamed[-CensoredInd] <- 1 

 return(EventVec_Renamed)
}
