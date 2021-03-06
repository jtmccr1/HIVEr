% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_munging.R
\name{getting_tp}
\alias{getting_tp}
\title{Finding transmission pairs}
\usage{
getting_tp(meta_one)
}
\arguments{
\item{meta_one}{A data frame containing one row for each ENROLLID in the data set.}
}
\value{
A data frame with one row for each possible transmission pair
}
\description{
Here we take a data frame and look for all cases where two individuals are
sick within a week of each other. This does not compute wether or not the
pairs match our criteria yet
}
\examples{
one_meta<-only_one(small_meta)
getting_tp(one_meta)
}
