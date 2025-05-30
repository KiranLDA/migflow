% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomSTAR.R
\name{randomSTAR}
\alias{randomSTAR}
\title{Generate a random network}
\usage{
randomSTAR(
  pop = 1e+05,
  toplot = TRUE,
  nbreeding = 8,
  nwintering = 4,
  nstop = 40,
  minforks = 2,
  maxforks = 5,
  anadromous = TRUE
)
}
\arguments{
\item{pop}{population size flowing through network}

\item{toplot}{TRUE/FALSE to determine whether the output is plotted or not}

\item{nbreeding}{number of breeding sites.}

\item{nwintering}{number of nonbreeding residency/wintering sites.}

\item{nstop}{number of sites used duting migratory stopover}

\item{minforks}{minimum number of forks in a river branch}

\item{maxforks}{maximum number of forks in a river branch}

\item{anadromous}{TRUE or FALSE. If TRUE, then population breeds at the tips of network branches, otheriwe, it breeding "at sea"}
}
\value{
a list containting the network which was randomly generated,
the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
}
\description{
This function generates a random migratory network
}
\examples{
par(mfrow=c(1,1))
randomSTAR(nbreeding = 6, toplot=TRUE)



}
