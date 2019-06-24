#' Generate a random network
#'
#' @description This function generates a random migratory network
#'
#' @param nsites number of sites to be generated in the network
#' @param nbreeding number of breeding sites. If all branch ends should be breeding sites then use "ALL"
#' @param toplot TRUE/FALSE to determine whether the output is plotted or not
#' @param pop population size flowing through network
#' @param minforks minimum number of forks in a river branch
#' @param maxforks maximum number of forks in a river branch
#' @param anadromous TRUE or FALSE. If TRUE, then population breeds at the tips of network branches, otheriwe, it breeding "at sea"
#'
#' @return a list containting the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#' par(mfrow=c(1,1))
#' randomRIVER(nsites = 20, toplot=TRUE)
#' randomRIVER(nsites = 50, toplot=TRUE, nbreeding="ALL")
#'
#'
#' @import igraph
#' @importFrom stats runif
#' @export
randomRIVER <- function(nsites = 50,
                        pop = 100000,
                        toplot= TRUE,
                        nbreeding = 8,
                        minforks = 2,
                        maxforks = 5,
                        anadromous = TRUE){

  # For testing purposes
  # nsites = 50
  # pop = 100000
  # toplot= TRUE
  # minforks = 2
  # maxforks = 5
  # anadromous = TRUE

  # Create a fake list of sites where animals were seen at, with latitude, longitude and number of anumals seen there
  site_list <- data.frame(Lat= runif(nsites, min=-20, max=40),
                          Lon= runif(nsites, min=-10, max=10),
                          Pop=runif(nsites, min=500, max=10000))

  #sort according to latitude
  site_list = site_list[order(site_list$Lat, decreasing=T),]
  site_list$Site= 1:nsites

  # add a dummy breeding and wintering site
  site_list<- rbind(
    c(41,0,pop,0),
    site_list,
    c(-21,0,pop,9999))

  # create a distance matrix based on these data
  dist <- point2DIST(site_list)
  rownames(dist)[1] <- "supersource"
  colnames(dist)[1] <- "supersource"
  rownames(dist)[nsites+2] <- "supersink"
  colnames(dist)[nsites+2] <- "supersink"

  # Make the network directed
  network <- directedNET(max(dist[2:(nsites+1), 2:(nsites+1)]) - dist[2:(nsites+1), 2:(nsites+1)], include_diagonal = TRUE)
  network[network==0] <- NA

  i=1
  forks <- round(runif(1,1,2))
  sorted <- network[i,]
  sorted[which(!is.na(sorted))] <- sort(network[i,which(!is.na(sorted))],index.return = TRUE)$ix
  keep <- which (sorted >= (max(sorted,na.rm=TRUE) - (forks-1)))
  # val= runif(length(keep),0,1)
  network[i,] <- NA
  network[i,keep] <- 1# val/(sum(val))

  i=2
  forks <- round(runif(1,minforks,maxforks))
  sorted <- network[i,]
  sorted[which(!is.na(sorted))] <- sort(network[i,which(!is.na(sorted))],index.return = TRUE)$ix
  keep <- which (sorted >= (max(sorted,na.rm=TRUE) - (forks-1)))
  already_inflowing <- which(!is.na(network[i-1,]))
  keep <- keep[which(!(keep %in% already_inflowing))]
  # val= runif(length(keep),0,1)
  network[i,] <- NA
  network[i,keep] <- 1# val/(sum(val))
  if (sum(network[,i], na.rm=TRUE) == 0){
    network[i-1,i]<- 1
  }


  for (i in 3:nrow(network)){
    forks <- round(runif(1,minforks,maxforks))
    sorted <- network[i,]
    sorted[which(!is.na(sorted))] <- sort(network[i,which(!is.na(sorted))],index.return = TRUE)$ix

    if(any(duplicated(sorted[!is.na(sorted)]))){
      idx <- which(any(duplicated(sorted[!is.na(sorted)])) == TRUE)
      sorted[!is.na(sorted)][idx:length(sorted)] <- sorted[!is.na(sorted)][idx:length(sorted)]+1
    }

    keep <- which (sorted >= (max(sorted,na.rm=TRUE) - (forks-1)))
    already_inflowing <- which(apply(network[1:i-1,],2,sum,na.rm=TRUE)>0)
    keep <- keep[which(!(keep %in% already_inflowing))]
    # val= runif(length(keep),0,1)
    network[i,] <- NA
    network[i,keep] <- 1# val/(sum(val))

    if (sum(network[,i], na.rm=TRUE) == 0){
      network[i-1,i]<- 1
    }
  }

  site_list$Site[site_list$Site == 0]<- "supersource"
  site_list$Site[site_list$Site == 9999]<- "supersink"

  #specify the sinks
  sinks = which(apply(network,1, function(x) sum(which(x>0)))==0)
  if(nbreeding != "ALL") sinks = sinks[which(sort(dist[sinks,"supersink"],index.return = TRUE)$ix <= nbreeding)]

  #Add supersource and sink nodes
  network <- addSUPERNODE(network, sources= site_list$Site[2],
                          sinks = sinks)

  network <-ifelse(network == Inf, pop, network)
  network[,"supersink"] <- ifelse(network[,"supersink"] == pop, max(dist[,"supersink"])-dist[,"supersink"], network[,"supersink"])
  network[,"supersink"] <- network[,"supersink"] /sum(network[,"supersink"] )#*pop

  network[is.na(network)] = 0
  network <- network * pop

  weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

  flow = max_flow(weight, source = V(weight)["supersource"],
                  target = V(weight)["supersink"], capacity = E(weight)$weight )

  nodes = get.edgelist(weight, names=TRUE)
  nodes = as.data.frame(nodes)
  nodes$flow = flow$flow
  nodes$Lat_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lat[site_list$Site %in% nodes[i,1]])))
  nodes$Lon_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lon[site_list$Site %in% nodes[i,1]])))
  nodes$Lat_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lat[site_list$Site %in% nodes[i,2]])))
  nodes$Lon_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lon[site_list$Site %in% nodes[i,2]])))


  neti <- weight
  E(neti)$weight <- flow$flow
  network <- as.matrix(as_adjacency_matrix(neti, attr="weight"))
  # neti[neti== "."] <- flow$flow

  if (toplot == TRUE){
    weight <- delete.vertices(weight, c(1, length(V(weight))))
    index = which(nodes$V1 != "supersource" & nodes$V2 != "supersink")
    nodeindex = which(nodes$V1 != "supersource")
    nodes[nodeindex,]
    sizes <- unlist(lapply(1:nsites ,function(x) sum(nodes$flow[nodeindex][nodes$V1[nodeindex]==x])))

    Bcols = c("black",ifelse(network[3:(nsites+1),"supersink"] > 0, "orange", "royalblue4"))
    lyt = layout_with_kk(weight)

    plot(weight, layout= lyt, edge.width = ((flow$flow[index]/pop)*20), edge.arrow.mode=0,
         edge.color = "royalblue4",
         vertex.color =  Bcols,
         vertex.label="",
         vertex.size = ((sizes/pop)*20)+5)
  }

  site_list$Pop <- c(pop,sizes, pop)
  site_list$Lon[2:(nrow(site_list)-1)] <- lyt[,1]
  site_list$Lat[2:(nrow(site_list)-1)] <- lyt[,2]

  if (anadromous == TRUE){
    colnames(network)[1]<- rownames(network)[1] <- "supersink"
    colnames(network)[ncol(network)]<- rownames(network)[nrow(network)] <- "supersource"
    mirror_network <- network # matrix(rev(network),ncol=52)
    mirror_network <- mirror_network[, ncol(mirror_network):1]
    mirror_network <- mirror_network[nrow(mirror_network):1, ]
    mirror_network <- t(mirror_network)
    SMnet <- mirror_network
    NMnet <- network
  } else {
    mirror_network <- network # matrix(rev(network),ncol=52)
    mirror_network <- mirror_network[, ncol(mirror_network):1]
    mirror_network <- mirror_network[nrow(mirror_network):1, ]
    mirror_network <- t(mirror_network)
    SMnet <- network
    NMnet <- mirror_network
  }


  #---------------------------------------------
  # Join the two networks
  #---------------------------------------------

  topright=matrix(0,nrow(SMnet),ncol(NMnet)+1)
  colnames(topright) <- c("NB",colnames(NMnet))

  bottomleft=matrix(0, nrow(NMnet)+1, ncol(SMnet))
  rownames(bottomleft) <- c("NB",rownames(NMnet))

  bottomright = cbind(rep_len(0,nrow(NMnet)),NMnet)
  bottomright = rbind(rep_len(0,ncol(bottomright)),bottomright)
  # bottomright[1,1] = 1

  top = cbind(SMnet,topright)
  bottom = cbind(bottomleft,bottomright)

  network=rbind(top,bottom)

  colnames(network)=c(paste0("S",colnames(SMnet)),"NB",paste0("N",colnames(NMnet)))
  rownames(network)=c(paste0("S",rownames(SMnet)),"NB",paste0("N",rownames(NMnet)))

  network["Ssupersink","NB"]<- pop
  network["NB","Nsupersink"]<- pop

  return(list( network = network,
               # tracks = tracks,
               sites = site_list))

}


