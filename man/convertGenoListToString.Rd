% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{convertGenoListToString}
\alias{convertGenoListToString}
\title{Convert from genotype list to string format}
\usage{
convertGenoListToString(genoList, genopheno)
}
\arguments{
\item{genoList}{A genotype representation in list format.}

\item{genopheno}{The genetic system that the genotype is part of}
}
\value{
Genotype in string format
}
\description{
Genotypes are internally represented as lists but specified as strings by the user.
This function converts from the former into the latter format.
}
\keyword{internal}
