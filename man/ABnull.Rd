% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABnull.R
\name{ABnull}
\alias{ABnull}
\title{Run model that considers no accumulation of epimutations (ABnull)}
\usage{
ABnull(pedigree.data, out.dir, out.name)
}
\arguments{
\item{pedigree.data}{Generation table name, you can find sample file in}

\item{out.dir}{outputdirectory}

\item{out.name}{name of file}
}
\value{
ABnull RData file.
}
\description{
Run model that considers no accumulation of epimutations (ABnull)
}
\examples{
#Get some toy data
inFile <- readRDS(system.file("extdata/dm/","output.rds", package="AlphaBeta"))
pedigree <- inFile$Pdata
out.name <- "CG_global_estimates_ABnull"
out <- ABnull(pedigree.data = pedigree,
                  out.dir=getwd(),
                  out.name=out.name)

summary(out)

}
