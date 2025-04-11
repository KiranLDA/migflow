#' Generate a random network
#'
#' @description This function generates a random migratory network
#'
#' @param nbreeding number of breeding sites.
#' @param nwintering number of nonbreeding residency/wintering sites.
#' @param nstop number of sites used duting migratory stopover
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
#' randomRIVER(nbreeding = 6, toplot=TRUE)
#'
#' @import igraph
#' @importFrom stats runif
#' @export
randomRIVER <- function(nwintering=4,
                        nbreeding=10,
                        nstop = 40,
                        pop = 100000,
                        toplot= TRUE,
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

  nsites <- nwintering + nbreeding + nstop

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
  if (anadromous==TRUE){
    # if(nbreeding != "ALL")
    sinks = sinks[which(sort(dist[sinks,"supersink"],index.return = TRUE)$ix <= nbreeding)]
    site_list$B <- 0
    site_list$B[sinks+1] <- 1
  }else{
    # if(nwintering != "ALL")
    sinks = sinks[which(sort(dist[sinks,"supersink"],index.return = TRUE)$ix <= nwintering)]
    site_list$NB <- 0
    site_list$NB[sinks+1] <- 1
  }


  # Add supersource and sink nodes
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

  nodes = as_edgelist(weight, names=TRUE)
  nodes = as.data.frame(nodes)
  nodes$flow = flow$flow
  nodes$Lat_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lat[site_list$Site %in% nodes[i,1]])))
  nodes$Lon_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lon[site_list$Site %in% nodes[i,1]])))
  nodes$Lat_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lat[site_list$Site %in% nodes[i,2]])))
  nodes$Lon_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(site_list$Lon[site_list$Site %in% nodes[i,2]])))


  # neti <- weight
  # E(neti)$weight <- flow$flow
  # network <- as.matrix(as_adjacency_matrix(neti, attr="weight"))
  # # neti[neti== "."] <- flow$flow


  weight <- delete.vertices(weight, c(1, length(V(weight))))
  index = which(nodes$V1 != "supersource" & nodes$V2 != "supersink")
  nodeindex = which(nodes$V1 != "supersource")
  nodes[nodeindex,]


  sizes <- apply(network,1,sum)

  if (anadromous==TRUE){
    index = as.numeric(names(sort(sizes[2:(length(sizes)-1)], decreasing =TRUE)))[1:nwintering]
    site_list$NB <- 0
    site_list$NB[index+1] <- 1
  } else{
    index = as.numeric(names(sort(sizes[2:(length(sizes)-1)], decreasing =TRUE)))[1:nbreeding]
    site_list$B <- 0
    site_list$B[index+1] <- 1
  }

  sizes <- unlist(lapply(1:nsites ,function(x) sum(nodes$flow[nodeindex][nodes$V1[nodeindex]==x])))



  site_list$NM <- site_list$SM <- 0
  site_list$NM[site_list$B == 0 & site_list$NB == 0] <- 1
  site_list$SM[site_list$B == 0 & site_list$NB == 0] <- 1

  # Bcols = c("black",ifelse(network[3:(nsites+1),"supersink"] > 0, "orange", "royalblue4"))
  lyt = layout_with_kk(weight)

  # if (toplot == TRUE){
  #     plot(weight, layout= lyt, edge.width = ((flow$flow[index]/pop)*20), edge.arrow.mode=0,
  #          edge.color = "royalblue4",
  #          vertex.color =  Bcols,
  #          vertex.label="",
  #          vertex.size = ((sizes/pop)*20)+5)
  #   }

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


  weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

  # run the population through the network a forst time
  flow = max_flow(weight, source = V(weight)["Ssupersource"],
                  target = V(weight)["Nsupersource"], capacity = E(weight)$weight)

  sites <-site_list
  if (toplot == TRUE){
    # plot flow network
    nodes = as_edgelist(weight, names=TRUE)
    nodes = as.data.frame(nodes)
    nodes$flow = flow$flow
    nodes$V1 <- substring(nodes$V1, 2)
    nodes$V2 <- substring(nodes$V2, 2)

    nodes = nodes[nodes$V1 != "supersource" & nodes$V1 != "supersink" & nodes$V2 != "supersource" & nodes$V2 != "supersink" ,]

    nodes$Lat_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lat[sites$Site %in% nodes[i,1]])))
    nodes$Lon_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lon[sites$Site %in% nodes[i,1]])))
    nodes$Lat_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lat[sites$Site %in% nodes[i,2]])))
    nodes$Lon_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lon[sites$Site %in% nodes[i,2]])))

    # library(shape)
    # par(mfrow=c(1,1))
    # par(mar=c(4,4,4,4))
    index=2:(nrow(sites)-1)
    if(toplot ==TRUE){
    plot(sites$Lon[index], sites$Lat[index], pch=16,
         cex=0, xlab="", ylab="", xaxt="n", yaxt = "n",
         frame.plot=FALSE)
      }

    index=1:nrow(nodes)
    segments(x0 = nodes$Lon_from[index],
             y0 = nodes$Lat_from[index],
             x1 = nodes$Lon_to[index],
             y1 = nodes$Lat_to[index],
             col= "black",
             lwd=(nodes$flow[index]/(max(nodes$flow)))*30)


    # sort sites by flow
    nodeflow = merge(aggregate(nodes$flow, by=list(Category=as.character(nodes$V1)), FUN=sum),
                     aggregate(nodes$flow, by=list(Category=as.character(nodes$V2)), FUN=sum), all=T)
    nodeflow$x = as.numeric(nodeflow$x)
    nodeflow = data.frame( unique(as.matrix(nodeflow[ , 1:2 ]) ))
    nodeflow$x = as.numeric(as.character(nodeflow$x))
    nodeflow = nodeflow[nodeflow$Category != "supersource" & nodeflow$Category != "supersink",]

    # make sure it is numeric
    nodeflow$Category = as.numeric(as.character(nodeflow$Category))

    # plot sites
    nodeflowplot = nodeflow[order(nodeflow$Category),]
    nodeflowplot$x[1] <-nodeflowplot$x[1]+nodeflowplot$x[1]
    # nodeflowplot = nodeflow[order(nodeflow$Category),]
    index=as.numeric(nodeflowplot$Category)+1
    colorz = ifelse(sites$B[index]==1,"royalblue",ifelse(sites$NB[index]==1,"orange","gray"))

    if(toplot ==TRUE){
    points(sites$Lon[index],
           sites$Lat[index],
           pch=21,
           cex=(((nodeflowplot$x)/
                   as.numeric(max(nodeflowplot$x)))+0.4)*4,
           bg=colorz , col="black")
    }
  }

  return(list( network = network,
               # tracks = tracks,
               sites = site_list))

}


