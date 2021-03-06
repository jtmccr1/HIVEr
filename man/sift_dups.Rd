% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing_snv.R
\name{sift_dups}
\alias{sift_dups}
\title{Sift duplicate calls
Take in a data frame of one mutation identified in 2 replicates and
choose the replicate with the better coverage.}
\usage{
sift_dups(df)
}
\arguments{
\item{df}{a data frame of nrow(2)}
}
\description{
Sift duplicate calls
Take in a data frame of one mutation identified in 2 replicates and
choose the replicate with the better coverage.
}
\examples{
sift_dups(subset(variants,mutation=="PB1_G530A"))
}
