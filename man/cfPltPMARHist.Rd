% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plottingFunctions.R
\name{cfPltPMARHist}
\alias{cfPltPMARHist}
\alias{cfPltFreqHist}
\title{Plot histograms of PMAR values by genotype}
\usage{
cfPltFreqHist(marVec, depVec, gtpVec, legend = FALSE, ff = NULL)
}
\arguments{
\item{marVec}{Integer, minor allele read counts}

\item{depVec}{Integer, total read depths}

\item{gtpVec}{Character, genotype pairs, e.g. "ABab"}

\item{legend}{logical of length 1, add a genotype legend when TRUE}

\item{ff}{Numeric of length 1, fetal fraction}
}
\description{
Plot histograms of PMAR values by genotype
}
\details{
When ff is not NULL, add normal approximation lines.
}
