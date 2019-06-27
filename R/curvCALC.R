#' Generate a random network
#'
#' @description This function generates a random migratory network
#'
#'
#' @param x value for x axis, typically number or proportion of sites lost
#' @param y value for y axis, typically population lost or proportion of population lost
#' @param digits number of digits to round the result to
#'
#' @return a list containting the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#'
#' pop=100000
#' rand_net = randomNET(nsites=100,nbreeding=3, nwintering=3,pop=pop,
#'                      mean_dist = 8000, sd_dist=9000, toplot=FALSE)
#'
#' # priotise sites according to flow through network
#' prioritisation = prioritiseFLOW(rand_net$network, rand_net$sites, toplot=FALSE, returnNET = FALSE)
#'
#' y=(prioritisation$prioritisation$Pop_Flow/pop)
#' x=((1:length(prioritisation$prioritisation$Pop_Flow))/
#'     length(prioritisation$prioritisation$Pop_Flow))
#'
#'
#' # calculate area under curve
#' AUC = curvCALC(x,y)
#' par(mfrow=c(1,1), mar=c(4,4,1,1))
#' plot(x,y,type="o",pch=20, main = AUC)
#'
#' @export
curvCALC <- function (x,y,digits=3){
  # library(ROCR)
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  return(round(abs(sum(dx * my)), digits = digits))
}
