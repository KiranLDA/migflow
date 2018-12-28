#' Only allow the network to flow from source to sink
#'
#' @description This function replaces the bottom side of the distance network with zeroes. This ensures that everything can only go in one direction. i.e. from code example A to B, A to C, A to D, B to C, B to D and C to D
#'
#' @param dist This is a distance matrix between two sites
#' @param include_diagonal `TRUE` means the diagonal is replaced (such that animals cannot remain in a site) and `FALSE` means they are not
#' @param replaceWith Replace lower part of distance matrix with e.g. zeros or NAs, default is set to 0
#'
#' @return A distance matrix in kilometres between all pairs sites
#'
#' @examples
#' dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
#' dist <- point2DIST(dta)
#' directedNET(dist)
#'
#'
#' @export
directedNET <- function(dist, include_diagonal= TRUE, replaceWith = 0){
  dist[lower.tri(dist, diag = include_diagonal)] <- replaceWith
  return(dist)
}
