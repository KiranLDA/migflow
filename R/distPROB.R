#' Calculate the distance probability
#'
#' @description This function calculates the probability of travelling from one point to the next given the distance between the points and the distance the species would normally travel in one go and
#'
#' @param tracks this is a series of distances travelled - for instance a distance that a bird has been observed to fly between two sites
#' @param dist This is a distance matrix between two sites
#' @param adjust this is a parameter for estimating the Kernal density of the tracked data. see '?density' for more details. The default is 2, but increasing will flatten and decreasing will add peaks.
#' @param plot Logical (true or false) for deciding whether or not to plot the kernal density distribution
#'
#' @return A matrix of probabilities showing how likely a bird can go from one site to the next given the distance between them
#'
#' @examples
#' tracks <- rnorm(10, 500, 200)
#' dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
#' dist <- point2DIST(dta)
#' distPROB(tracks, dist, adjust=2, plot=TRUE)
#'
#' @importFrom stats approxfun density
#' @export
distPROB <- function(tracks, dist, adjust = 2, plot= TRUE){
  prob=approxfun(density(tracks,adjust=adjust))
  if (plot) plot(density(tracks,adjust=adjust), xlab= "Distance")
  dist_prob = apply(dist,c(1,2),function(x) prob(x))
  dist_prob[is.na(dist_prob)]=0
  return(dist_prob)
}
