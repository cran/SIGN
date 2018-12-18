#' ExpPheno_Categorize is a function for grouping samples 
#' based on their survival to 3 groups of poor, good, and intermediate 
#' @param ExpMeta_List List containing expression matrix and metadata matrix
#' @param Time_ID Index of time to death in metadata matrix
#' @param Event_ID Index of event in metadata matrix
#' @param Mad_Factor Threshold of mad in time to death values to determine poor survival group
#' @param MinNum_ExClass Minimum number of samples that has to be kept in poor and good group 
#' (if number of samples is lower than this threhold, more samples will be addedd in order of survival)
#' @param Expression_Log2 Parameter for gene expsression value transformation to logarithmic scale (log2(expression value+1))
#' @return List of expression matrices, and time to event as well as event for the patients
#'  within each category of poor, intermediate or good survival
#' @importFrom stats cov median
#' @importFrom stats cov mad
#' @importFrom stats cov na.omit

ExpPheno_Categorize <- function(ExpMeta_List, Time_ID, Event_ID, Mad_Factor, MinNum_ExClass, Expression_Log2=FALSE){

 ExpMat <- ExpMeta_List[["Expression"]]
 MetaMat <- ExpMeta_List[["Meta"]]
 
 TimeVec <- as.numeric(as.character(MetaMat[,Time_ID]))
 EventVec <- as.character(MetaMat[,Event_ID])
#################
 LowSurvInd <- which(log10(TimeVec) < (median(log10(TimeVec))- Mad_Factor*mad(log10(TimeVec))) & (EventVec == "1" | tolower(EventVec) == "deceased"))

 if(length(LowSurvInd) < MinNum_ExClass){
  LowSurvInd <- sort(TimeVec, decreasing = F, index.return = T)[[2]][c(1:MinNum_ExClass)]
 }
#################
 HighSurvInd <- which(TimeVec > (median(TimeVec) + mad(TimeVec)))
 if(length(HighSurvInd) < MinNum_ExClass){
  HighSurvInd <- sort(TimeVec, decreasing = T, index.return = T)[[2]][c(1:MinNum_ExClass)]
 }
#################
 OverlapInd <- intersect(LowSurvInd, HighSurvInd)
 if(length(OverlapInd) > 0){
  stop("There is not enough sample to separte poor and good classes!")
 }
#################
 MedSurvInd <- seq(1,length(as.numeric(TimeVec)))[-c(HighSurvInd, LowSurvInd)]
################
 MaxExp <- max(na.omit(as.numeric(ExpMat)))
 if(Expression_Log2){
  if(MaxExp > 100){
   ExpMat <- log2(ExpMat+1)
  }
 }else{
  if(MaxExp < 100){
   ExpMat <- (2^ExpMat-1)
  }
 }
#################
 LowSuvMat <- ExpMat[,LowSurvInd]
 MedSuvMat <- ExpMat[,MedSurvInd]
 HighSuvMat <- ExpMat[,HighSurvInd]

 ExpList_Categorized <- list(LowSuvMat, MedSuvMat, HighSuvMat)
 names(ExpList_Categorized) <- c("poor", "intermediate", "good")

 SurvivalDayMat <- list(TimeVec[LowSurvInd], TimeVec[MedSurvInd], TimeVec[HighSurvInd])
 names(SurvivalDayMat) <- c("poor", "intermediate", "good")
 
 VitalStatusMat <- list(EventVec[LowSurvInd], EventVec[MedSurvInd], EventVec[HighSurvInd])
 names(VitalStatusMat) <- c("poor", "intermediate", "good")

 ExpPheno_Categorized <- list(ExpList_Categorized, SurvivalDayMat, VitalStatusMat)
 names(ExpPheno_Categorized) <- c("ExpList", "TimeList", "EventList")

 return(ExpPheno_Categorized)
}
