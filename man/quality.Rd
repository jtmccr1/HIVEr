% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing_snv.R
\name{quality}
\alias{quality}
\title{Identify quality variants}
\usage{
quality(df)
}
\arguments{
\item{df}{a data frame containing all base calls associated with a specimen or isolate.}
}
\description{
This function takes in all variant calls associated with a specimen.
If the sample was sequenced in duplicate that must be apparent by a 'run' column
that specifies the unique names of the run columns. Samples with gc_ul above 1e5 don't
need to be run in duplicate. If the sample has a lower titer but was not sequenced twice this
function will through an error.
}
\examples{
quality(variants)

}
