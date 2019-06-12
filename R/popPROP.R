#' Calculate the proportion of population that can flow into each site
#'
#' @description This function calculates determines the proportion of the population flowing into a site - this is used to later parameterise preference
#'
#' @param network network of connections between sites
#' @param population This is the total size of the populaiton or number of individuals allowed to flow through the network
#'
#' @return A matrix containing the proportion of the population able to flow into that
#'
#' @examples
#' dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
#' network <- point2DIST(dta)
#' popPROP(network, 300)
#'
#'
#'
#' @export
popPROP <- function(network, population){
  prop = t(apply(network,1,
                   function(x) x[which(!is.na(x))]/
                     sum(x,na.rm=TRUE)))
  prop[1,] = population*prop[1,]
  # output = prop
  # output[] = c(prop[1,],unlist(lapply(2:nrow(prop), function(i) {prop[i,] <- sum(prop[1:(i-1),i])*prop[i,]})))
  for(i in 2:nrow(prop)){
    prop[i,] = sum(prop[1:(i-1),i])*prop[i,]
  }
  prop[is.na(prop)] <- 0
  return(prop)
}

