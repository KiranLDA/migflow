#' Generate a random network
#'
#' @description This function generates a random migratory network
#'
#'
#' @param nsites number of sites to be generated in the network
#' @param nedges number of edges in the network. If all sites should be connected then use nedges="ALL"
#' @param pop population size flowing through network
#' @param mean_dist average distance the species can travel (this can be estimated from real data)
#' @param sd_dist standard deviation of the distances a species can travel (this can be estimated from real data)
#' @param Latrange geographic range of species, by default latitude restricted to c(-20,40),
#' @param Lonrange geographic range of species, by default longitude restricted to c(-10,10),
#' @param Poprange min and max population sizes to be randomly generated, by default c(100,10000)
#' @param plot TRUE/FALSE to determine whether the output is plotted or not
#'
#' @return a list containting the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#'
#' pop=100000
#' rand_net = randomNET(nsites=100,nbreeding=3, nwintering=3,pop=pop, mean_dist = 8000, sd_dist=9000, plot=F)
#'
#' # priotise sites according to flow through network
#' prioritisation = prioritiseFLOW(rand_net$network, rand_net$sites, plot=F)
#'
#' y=(prioritisation$prioritisation$Pop_Flow/pop)
#' x=(1:length(prioritisation$prioritisation$Pop_Flow))/length(prioritisation$prioritisation$Pop_Flow)
#'
#'
#' # calculate area under curve
#' AUC = curvCALC(x,y)
#' par(mfrow=c(1,1), mar=c(4,4,1,1))
#' plot(x,y,type="o",pch=20, main = AUC)
#'
#' @importFrom base diff length sum abs round
#' @export
curvCALC <- function (x,y,digits=3){
  # library(ROCR)
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  return(round(abs(sum(dx * my)), digits = digits))
}
