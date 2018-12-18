#' Pathway_Grouping is a function to make a pathway list from files containing genes within each pathway
#' @param PathwayDir Path of directory including the files of pathways
#' @param Pattern Pattern should be used to select the files of pathway genes from PathwayDir
#' @return List of genes within the pathway

Pathway_Grouping <- function(PathwayDir, Pattern){
 
 PathwayFiles <- list.files(PathwayDir, Pattern)

 if(grepl("rds", Pattern)){
  PathwayList <- list()
  for(FileIter in 1:length(PathwayFiles)){
   PathwayList[[FileIter]] <- readRDS(paste(PathwayDir, PathwayFiles[FileIter],sep = "", collapse = ""))
  }
  names(PathwayList) <- gsub(Pattern, "", PathwayFiles)
 }else{
  stop("Please use Pathways_Preprocess_MSigDB.R script and save the files as rds!")
 }

 return(PathwayList)
}
