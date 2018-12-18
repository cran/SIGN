#' Survival_Stats is a function for building cox model using all the features and each feature as a separate model 
#' @param ScoreMat Matrix of feature values used for survival predition
#' @param TimeVec Vectore of time to death of samples (patients)
#' @param EventVec Vector of events for the samples (patients) as being dead or censored
#' @return A list containing summary of a cox model using all of the features and separate cox models for each feature
#' @importFrom survival Surv
#' @importFrom survival coxph


Survival_Stats <- function(ScoreMat, TimeVec, EventVec){

 TimeVec <- as.numeric(TimeVec)
 EventVec <- as.numeric(EventVec)

 
 ScoreMat[,which(grepl("wilcox", tolower(colnames(ScoreMat))))] <- -log10(ScoreMat[,which(grepl("wilcox", tolower(colnames(ScoreMat))))])
 ScoreMat[which(ScoreMat == Inf)] <- 10*max(ScoreMat[which(ScoreMat != Inf)]) 

 Signature <- as.data.frame(princomp(ScoreMat)$scores[,c(1:ceiling(ncol(ScoreMat)/2))])

 SurvObject <- Surv(as.numeric(TimeVec), as.numeric(EventVec))

 Signature$SurvObj <- SurvObject

 Cox_Summary <- summary(coxph(formula = SurvObj ~ ., data = Signature))
 
 CoxFeatures <- list()
 for(FeatIter in 1:ncol(ScoreMat)){
  CoxFeatures[[FeatIter]]  <- summary(coxph(formula = SurvObject ~ as.numeric(ScoreMat[,FeatIter])))
 }

 names(CoxFeatures) <- colnames(ScoreMat)

 Stat_List <- list(Cox_Summary, CoxFeatures)
 names(Stat_List) <- c("Cox_All", "Cox_SepFeatures")

 return(Stat_List)
}
