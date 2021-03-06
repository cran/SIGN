% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIGN_Ensemble_SimCal.R
\name{SIGN_Ensemble_SimCal}
\alias{SIGN_Ensemble_SimCal}
\title{SIGN_Ensemble_SimCal is a function for Generating list fo similarities based on 
different pathway quantification methods and similarity measures}
\usage{
SIGN_Ensemble_SimCal(ExpList, RefClassID, TestClassID, GeneID, PathwaySets)
}
\arguments{
\item{ExpList}{List of expression matrices for different groups of samples
used in the centroid classification scheme}

\item{RefClassID}{Names of the matrices in the ExpList}

\item{TestClassID}{ID of a matrix in ExpList to be used as test set}

\item{GeneID}{Parameter to determine if gene annotations are provided as Symbols or EntrezIDs}

\item{PathwaySets}{List of pathways containing gene annotations for each pathways}
}
\value{
List of similarities identified in both gene and pathway level
}
\description{
SIGN_Ensemble_SimCal is a function for Generating list fo similarities based on 
different pathway quantification methods and similarity measures
}
