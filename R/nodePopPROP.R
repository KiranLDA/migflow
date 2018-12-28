#' Calculate the proportion of population that can flow into each site
#'
#' @description This function calculates determines the proportion of the population flowing into a site - this is used to later parameterise preference
#'
#' @param dta This is a dataframe. Columns of interest must be named: Site, Lat (for latitude), Lon (for longitude).
#' @param population This is the total size of the populaiton or number of individuals allowed to flow through the network
#'
#' @return A matrix containing the proportion of the population able to flow into that
#'
#' @examples
#' dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
#' nodePopPROP(dta, 300)
#'
#'
#'
#' @export
nodePopPROP <- function(dta, population = 1000){
  #create dataset to fill
  dist = point2DIST(dta)
  #get population information from each site
  for (i in 1:(dim(dist)[2])){
    dist[i,1:(dim(dist)[2])]= dta$Pop/population
  }
  dist[is.na(dist)] = 0
  return(dist)
}
