#' Generate a random network
#'
#' @description This function bridges one wintering node to another (rather than going through a supersink)
#'
#' @param network the network to be shortened
#' @param from the site for the network to be shortened from
#' @param to the site for the network to be shortened to
#'
#' @return a list containting the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#' par(mfrow=c(1,1))
#' net <- randomEIGHT(toplot=FALSE)$network
#' net <- shortenNET(net, from = "Ssupersink", to = "Nsupersink")
#'
#' @export
shortenNET <- function(network, from = "Ssupersink", to = "Nsupersink"){
  arrivalsites <- network[,from][network[,from]>0]
  for (i in 1:length(arrivalsites)){
    coli <- which(round(as.numeric(network[to,])) == round(as.numeric(arrivalsites[i])))
    rowi <- which(round(as.numeric(network[,from])) == round(as.numeric(arrivalsites[i])))
    network[rowi, coli] <- arrivalsites[i]
  }

  row_remove <- -(which(rownames(network) == from): which(rownames(network) == to))
  col_remove <- -(which(colnames(network) == from): which(colnames(network) == to))

  network <- network[row_remove,col_remove]
  return(network)
}
