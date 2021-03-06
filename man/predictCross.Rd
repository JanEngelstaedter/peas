% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossing.R
\name{predictCross}
\alias{predictCross}
\title{Predicts offspring distributions resulting from a cross between two genotypes.}
\usage{
predictCross(genopheno, mom, dad, equivalent = "phase",
  output = "both")
}
\arguments{
\item{genopheno}{genopheno object specifying the genetic setup.}

\item{mom}{Genotype of the mother.}

\item{dad}{Genotype of the father.}

\item{equivalent}{Argument determining which genotypes are to be treated as equivalent.
Possible options are "phase" (genotypes that differ only by phase are equivalent,
i.e. the order of alleles within a locus does not matter),
"origin" (phase is important but not the origin of alleles),
and "none" (no two genotypes are equivalent).}

\item{output}{Specifies what the function should return: 'genotypes' returns the distribution of
offspring genotypes, 'phenotypes' returns the distribution of offspring phenotypes, and
'both' returns both distributions.}
}
\value{
Depending on the output argument, a dataframe containing the distribution of
offspring genotypes, a dataframe containing the distribution of offspring phenotypes,
or both (combined in a list).
}
\description{
Predicts offspring distributions resulting from a cross between two genotypes.
}
\examples{
MendelsPeas <- newGenopheno(alleleNames = list(c('Y', 'y')))
MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'Y_', 'yellow')
MendelsPeas <- setPhenotypes(MendelsPeas, 'colour', 'yy', 'green')
predictCross(MendelsPeas, 'YY', 'yy')
predictCross(MendelsPeas, 'Yy', 'Yy')

}
