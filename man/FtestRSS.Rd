% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FtestRSS.R
\name{FtestRSS}
\alias{FtestRSS}
\title{Comparison of different models and selection of best model}
\usage{
FtestRSS(pedigree.select, pedigree.null)
}
\arguments{
\item{pedigree.select}{pedigree model.}

\item{pedigree.null}{ABnull pedigree.}
}
\value{
result of Ftest.
}
\description{
Comparison of different models and selection of best model
}
\examples{
## Get some toy data
file1 <- system.file("extdata/models/","ABneutral_CG_global_estimates.Rdata", package="AlphaBeta")
file2 <- system.file("extdata/models/","ABnull_CG_global_estimates.Rdata", package="AlphaBeta")
out <- FtestRSS(pedigree.select=file1,
                pedigree.null=file2)
}
