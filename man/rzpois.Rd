% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transmission_models.R
\name{rzpois}
\alias{rzpois}
\title{Randomly sample a zero truncated Poisson}
\usage{
rzpois(n, lambda)
}
\arguments{
\item{n}{number samples}

\item{lambda}{The mean of the untruncated Poisson}
}
\value{
a random sample
}
\description{
Randomly sample a zero truncated Poisson distribution given
the lambda of the underlying Poisson (that is of course truncated).
}
\examples{
rzpois(10,2.3)
}
