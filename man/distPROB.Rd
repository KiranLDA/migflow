% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distPROB.R
\name{distPROB}
\alias{distPROB}
\title{Calculate the distance probability}
\usage{
distPROB(tracks, dist, adjust = 2, plot = TRUE)
}
\arguments{
\item{tracks}{this is a series of distances travelled - for instance a distance that a bird has been observed to fly between two sites}

\item{dist}{This is a distance matrix between two sites}

\item{adjust}{this is a parameter for estimating the Kernal density of the tracked data. see '?density' for more details. The default is 2, but increasing will flatten and decreasing will add peaks.}

\item{plot}{Logical (true or false) for deciding whether or not to plot the kernal density distribution}
}
\value{
A matrix of probabilities showing how likely a bird can go from one site to the next given the distance between them
}
\description{
This function calculates the probability of travelling from one point to the next given the distance between the points and the distance the species would normally travel in one go and
}
\examples{
tracks <- rnorm(10, 500, 200)
dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
dist <- point2DIST(dta)
distPROB(tracks, dist, adjust=2, plot=TRUE)

}
