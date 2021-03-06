% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutational_model.R
\name{p0}
\alias{p0}
\title{Mutational model- non polmorphic}
\usage{
p0(mu, Ne, t)
}
\arguments{
\item{mu}{The mutation rate per nucleotide per cellular infectious cycle}

\item{Ne}{The effective population size}

\item{t}{Generations}
}
\description{
What's the probability the site still is at frequency 0?
This assumes anything below 0.001 (0.1%) is at 0. It yields 1- integral g_ft from 0.001 to 1
}
