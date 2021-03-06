% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Genes_SimCal.R
\name{Genes_SimCal}
\alias{Genes_SimCal}
\title{Genes_SimCal is a function to calculate similarity between a set of samples
 and 2 reference groups of samples}
\usage{
Genes_SimCal(ExpMat_Test, ExpMat_Ref1, ExpMat_Ref2, RefIDs, TestClassIter,
  SampleIter)
}
\arguments{
\item{ExpMat_Test}{Expression matrix for the test samples for which SIGN will 
indetify the similarity with the 2 reference sataset}

\item{ExpMat_Ref1}{Expression matrix for the 1st reference set fo samples}

\item{ExpMat_Ref2}{Expression matrix for the 2nd reference set fo samples}

\item{RefIDs}{Annotations corresponding to the 2 expression matrices 
(1st and 2nd names are associated with the 1st and 2nd expression matrix and )}

\item{TestClassIter}{Index to be matched with RefIDs for removal of test samples
from reference expression matrices}

\item{SampleIter}{Index of samples in the test expression matrix exist in 
referencece expression matrix 1 or 2}
}
\value{
Vector of similarity between the target samples and the 2 reference sets
}
\description{
Genes_SimCal is a function to calculate similarity between a set of samples
 and 2 reference groups of samples
}
