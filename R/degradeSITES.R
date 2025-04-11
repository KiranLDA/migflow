#' Gradually degrade sites
#'
#' @description Degrate sites gradually
#'
#' @param network adjacency network describing where birds are going
#' @param sites network sites including latitude and longitude, and population count at each site
#' @param degrade sites to degrade.Typically `c( "B", "STP", "NB")`. Can be B for breeding, STP for stopover sites, and/or NB for non-breeding residency
#' @param mindeg minimum percentage to degrade the site by (i.e. reduce carrying capacity by 1 percent)
#' @param maxdeg minimum percentage to degrade the site by (i.e. reduce carrying capacity by 10 percent)
#' @param iter number of iterations, if wanting the model to run until no population is left, then iter=NA
#' @param returnNET TRUE/FALSE whether or not the network is a return network or not
#' @param toplot TRUE/FALSE to determine whether the output is plotted or not
#' @param pop number of individuals in the population
#'
#' @return a list containting the prioritisation (as a list of the sites in the order they are removed), the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#' # with a return network
#' net <- randomPARALLEL(nwintering=5,mean_dist=0, sd_dist=1500, nbreeding=5)
#' network <- shortenNET(net$network, from = "Ssupersink", to = "Nsupersink")
#'
#' colnames(network)[1]<- rownames(network)[1] <- "Ssupersource"
#' colnames(network)[length(colnames(network))] <-  "Nsupersink"
#' rownames(network)[length(rownames(network))] <- "Nsupersink"
#' sites <- net$sites
#' degradeSITES(network,sites)
#'
#' @import igraph
#' @importFrom graphics par segments points
#' @importFrom stats aggregate
#' @importFrom grDevices colorRampPalette
#' @export
degradeSITES <- function(network,
                         sites,
                         degrade = c( "B", "STP", "NB"),
                         mindeg = 1,
                         maxdeg = 10,
                         toplot = TRUE,
                         pop = 100000,
                         iter=NA,
                         returnNET = TRUE){


  # for testing
  # net <- randomBOTTLE(nwintering=5,mean_dist=0, sd_dist=1500, nbreeding=5)
  # network <- shortenNET(net$network, from = "Ssupersink", to = "Nsupersink")
  #
  # colnames(network)[1]<- rownames(network)[1] <- "Ssupersource"
  # colnames(network)[length(colnames(network))] <-  "Nsupersink"
  # rownames(network)[length(rownames(network))] <- "Nsupersink"
  # sites <- net$sites



  if(returnNET == TRUE){

    #created a weigted igraph network
    weight <- graph_from_adjacency_matrix(network, mode="directed", weighted = TRUE)

    sourcename = rownames(network)[1]
    sinkname =  colnames(network)[length(colnames(network))]

    # run the population through the network a forst time
    flow = max_flow(weight, source = V(weight)[sourcename],
                    target = V(weight)[sinkname],
                    capacity = E(weight)$weight )


    nodes = as_edgelist(weight, names=TRUE)
    nodes = as.data.frame(nodes)
    nodes$flow = flow$flow #E(weight)$weight#
    # nodes3 <- nodes
    nodes$V1 <- substring(nodes$V1, 2)
    nodes$V2 <- substring(nodes$V2, 2)

    nodes$Lat_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lat[sites$Site %in% nodes[i,1]])))
    nodes$Lon_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lon[sites$Site %in% nodes[i,1]])))
    nodes$Lat_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lat[sites$Site %in% nodes[i,2]])))
    nodes$Lon_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lon[sites$Site %in% nodes[i,2]])))

    nodes2 = nodes[nodes$V1 != "supersource" & nodes$V2 != "supersink" ,]


    # sort sites by flow
    nodeflow = merge(aggregate(nodes$flow, by=list(Category=as.character(nodes$V1)), FUN=sum),
                     aggregate(nodes$flow, by=list(Category=as.character(nodes$V2)), FUN=sum), all=T)

    nodeflow$x = as.numeric(nodeflow$x)
    nodeflow = data.frame( unique(as.matrix(nodeflow[ , 1:2 ]) ))
    nodeflow$x = as.numeric(as.character(nodeflow$x))
    nodeflow = nodeflow[nodeflow$Category != "supersource" & nodeflow$Category != "supersink",]

    # specify the site to remove, in this case, it is the site which contributes less to population flow
    to_remove = which(nodeflow$x %in% min(nodeflow$x))[1]

    # calculate betweeness
    betweenness <- betweenness(weight, weights = flow$flow+0.0001)
    names(betweenness) <- substring(names(betweenness), 2)
    between <- c()
    for (n in as.character(nodeflow$Category[to_remove])){
      between <- c(between, sum(betweenness[which(names(betweenness) == n)]))
    }


    # empty dataset to store output
    prioritisation <- data.frame(Site=as.character(nodeflow$Category[to_remove]),
                                 Pop_Flow =  flow$value,#as.numeric(sites$Pop[1]),#
                                 Site_Flow =   nodeflow$x[to_remove],
                                 betweenness = between,
                                 average_path = average.path.length(weight))
    # make sure it is numeric
    nodeflow$Category = as.numeric(as.character(nodeflow$Category))

    if (toplot == TRUE){
      par(mar=c(2,2,2,2), mfrow=c(2,2))

      #plot empty points
      index=2:(nrow(sites)-1)
      plot(sites$Lon[index], sites$Lat[index], pch=16,
           cex=0, xlab="", ylab="", xaxt="n", yaxt = "n",
           frame.plot=FALSE)

      #plot edges
      index=1:nrow(nodes2)
      segments(x0 = nodes2$Lon_from[index],
               y0 = nodes2$Lat_from[index],
               x1 = nodes2$Lon_to[index],
               y1 = nodes2$Lat_to[index],
               col= "black",
               lwd=(nodes2$flow[index]/(max(nodes2$flow)))*15)

      #plot nodes
      index=2:(nrow(sites)-1)
      colorz = ifelse(sites$B[index]==1,"royalblue",ifelse(sites$NB[index]==1,"orange","gray"))
      points(sites$Lon[index],
             sites$Lat[index],
             pch=21,
             cex=((as.numeric(sites$Pop[index])/max(as.numeric(sites$Pop[index])))+0.4)*4,
             bg=colorz,
             col="black")
    }


    # store full network for results, before things get removed in the optimisation
    colnames(network) <- substring(colnames(network), 2)
    rownames(network) <- substring(rownames(network), 2)

    full_network = network
    full_site_list = sites

    #reduce the capacities of different types of nodes

    if ("B" %in% degrade){
      network[1,which(network[1,] >0)] <- unlist(lapply(which(network[1,] >0),
                                                        function(x) return(network[1,x] -(network[1,x] * runif(1,mindeg,maxdeg)/100))))
      network[,ncol(network)] <- rev(network[1,])
    }

    if("NB"  %in% degrade){
      rid <- which(rownames(network) %in% sites$Site[sites$NB == 1])[1:length(which(sites$NB == 1))]
      cid <- which(colnames(network) %in% sites$Site[sites$NB == 1])[(length(which(sites$NB == 1))+1):(length(which(sites$NB == 1))*2)]
      replace <- network[rid, cid]
      replace[replace>0] <- unlist(lapply(which(replace >0),
                                          function(x) return(replace[x] -(replace[x] * runif(1,mindeg,maxdeg)/100))))

      network[rid, cid] <- replace
      rm(rid)
      rm(cid)
      rm(replace)

    }
    if ("STP" %in% degrade){
      cid <- which(colnames(network) %in% which(sites$NM == 1 & sites$SM == 1))
      network[,cid]<-apply(network[,cid],2, function(x) x -(x * runif(1,mindeg,maxdeg)/100))
      rm(cid)

    }

    network[network<0]<- 0
    # prioritisation

    # net_remove = which(colnames(network) %in% as.character(nodeflow$Category[to_remove]))
    # network = network[-net_remove,-net_remove ]
    # sites = sites[!(sites$Site %in% as.character(nodeflow$Category[to_remove])),]
    #
    #
    if(is.na(iter)){
      while(flow$value > ceiling(pop*0.001)){

        weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

        sourcename = rownames(network)[1]
        sinkname =  colnames(network)[length(colnames(network))]

        # run the population through the network a forst time
        flow = max_flow(weight, source = V(weight)[sourcename],
                        target = V(weight)[sinkname],
                        E(weight)$weight)

        nodes = as_edgelist(weight, names=TRUE)
        nodes = as.data.frame(nodes)
        nodes$flow = flow$flow
        # nodes$V1 <- substring(nodes$V1, 2)
        # nodes$V2 <- substring(nodes$V2, 2)

        nodeflow = merge(aggregate(nodes$flow, by=list(Category=as.character(nodes$V1)), FUN=sum),
                         aggregate(nodes$flow, by=list(Category=as.character(nodes$V2)), FUN=sum), all=T)
        nodeflow$x = as.numeric(nodeflow$x)

        nodeflow = data.frame( unique(as.matrix(nodeflow[ , 1:2 ]) ))
        nodeflow$x = as.numeric(as.character(nodeflow$x))

        nodeflow = nodeflow[nodeflow$Category != "supersource" & nodeflow$Category != "supersink",]
        to_remove = which(nodeflow$x %in% min(nodeflow$x))[1]

        betweenness <- betweenness(weight, weights = flow$flow+0.0001)
        # names(betweenness) <- substring(names(betweenness), 2)
        between <- c()
        for (n in as.character(nodeflow$Category[to_remove])){
          between <- c(between, sum(betweenness[which(names(betweenness) == n)]))
        }



        prioritisation <- rbind(prioritisation,
                                data.frame(Site=as.character(nodeflow$Category[to_remove]),
                                           Pop_Flow =   flow$value,
                                           Site_Flow =   nodeflow$x[to_remove],
                                           betweenness = between,
                                           average_path = average.path.length(weight)))

        if ("B" %in% degrade){
          network[1,which(network[1,] >0)] <- unlist(lapply(which(network[1,] >0),
                                                            function(x) return(network[1,x] -(network[1,x] * runif(1,mindeg,maxdeg)/100))))
          network[,ncol(network)] <- rev(network[1,])
        }

        if("NB"  %in% degrade){
          rid <- which(rownames(network) %in% sites$Site[sites$NB == 1])[1:length(which(sites$NB == 1))]
          cid <- which(colnames(network) %in% sites$Site[sites$NB == 1])[(length(which(sites$NB == 1))+1):(length(which(sites$NB == 1))*2)]
          replace <- network[rid, cid]
          replace[replace>0] <- unlist(lapply(which(replace >0),
                                              function(x) return(replace[x] -(replace[x] * runif(1,mindeg,maxdeg)/100))))

          network[rid, cid] <- replace
          rm(rid)
          rm(cid)
          rm(replace)

        }
        if ("STP" %in% degrade){
          cid <- which(colnames(network) %in% which(sites$NM == 1 & sites$SM == 1))
          network[,cid]<-apply(network[,cid],2, function(x) x -(x * runif(1,mindeg,maxdeg)/100))
          rm(cid)

        }
        network[network<0]<- 0

      }

    } else{
      for(z in 1:iter){
        # if(flow$value > ceiling(pop*0.001)){

        weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

        sourcename = rownames(network)[1]
        sinkname =  colnames(network)[length(colnames(network))]

        # run the population through the network a forst time
        flow = max_flow(weight, source = V(weight)[sourcename],
                        target = V(weight)[sinkname],
                        E(weight)$weight)

        nodes = as_edgelist(weight, names=TRUE)
        nodes = as.data.frame(nodes)
        nodes$flow = flow$flow
        # nodes$V1 <- substring(nodes$V1, 2)
        # nodes$V2 <- substring(nodes$V2, 2)

        nodeflow = merge(aggregate(nodes$flow, by=list(Category=as.character(nodes$V1)), FUN=sum),
                         aggregate(nodes$flow, by=list(Category=as.character(nodes$V2)), FUN=sum), all=T)
        nodeflow$x = as.numeric(nodeflow$x)

        nodeflow = data.frame( unique(as.matrix(nodeflow[ , 1:2 ]) ))
        nodeflow$x = as.numeric(as.character(nodeflow$x))

        nodeflow = nodeflow[nodeflow$Category != "supersource" & nodeflow$Category != "supersink",]
        to_remove = which(nodeflow$x %in% min(nodeflow$x))[1]

        betweenness <- betweenness(weight, weights = flow$flow+0.0001)
        # names(betweenness) <- substring(names(betweenness), 2)
        between <- c()
        for (n in as.character(nodeflow$Category[to_remove])){
          between <- c(between, sum(betweenness[which(names(betweenness) == n)]))
        }



        prioritisation <- rbind(prioritisation,
                                data.frame(Site=as.character(nodeflow$Category[to_remove]),
                                           Pop_Flow =   flow$value,
                                           Site_Flow =   nodeflow$x[to_remove],
                                           betweenness = between,
                                           average_path = average.path.length(weight)))

        if ("B" %in% degrade){
          network[1,which(network[1,] >0)] <- unlist(lapply(which(network[1,] >0),
                                                            function(x) return(network[1,x] -(network[1,x] * runif(1,mindeg,maxdeg)/100))))
          network[,ncol(network)] <- rev(network[1,])
        }

        if("NB"  %in% degrade){
          rid <- which(rownames(network) %in% sites$Site[sites$NB == 1])[1:length(which(sites$NB == 1))]
          cid <- which(colnames(network) %in% sites$Site[sites$NB == 1])[(length(which(sites$NB == 1))+1):(length(which(sites$NB == 1))*2)]
          replace <- network[rid, cid]
          replace[replace>0] <- unlist(lapply(which(replace >0),
                                              function(x) return(replace[x] -(replace[x] * runif(1,mindeg,maxdeg)/100))))

          network[rid, cid] <- replace
          rm(rid)
          rm(cid)
          rm(replace)

        }
        if ("STP" %in% degrade){
          cid <- which(colnames(network) %in% which(sites$NM == 1 & sites$SM == 1))
          network[,cid]<-apply(network[,cid],2, function(x) x -(x * runif(1,mindeg,maxdeg)/100))
          rm(cid)

        }
        network[network<0]<- 0

        if (flow$value <= ceiling(pop*0.001)) break

      }
    }#}




    y = prioritisation$Pop_Flow/prioritisation$Pop_Flow[1]
    x = (1:length(prioritisation$Pop_Flow))/length(prioritisation$Pop_Flow)
    AUC = curvCALC(x,y)
    if (toplot == TRUE){
      par(mar=c(4,4,1,1))
      plot(prioritisation$Pop_Flow, type="o",pch=16,
           col = colorRampPalette(c("royalblue", "orange"))(nrow(full_site_list)-1)[1:(nrow(full_site_list)-1)],
           xlab="# Sites Removed",
           ylab="Population size",
           main=AUC)

      par(mar=c(4,4,1,1))
      plot(prioritisation$Site_Flow, c(0,abs(diff(prioritisation$Pop_Flow))), pch=16,
           col = colorRampPalette(c("royalblue", "orange"))(nrow(full_site_list)-1)[1:(nrow(full_site_list)-1)],
           xlab="Site carrying capacity",
           ylab="Change in population size when site lost")

      par(mar=c(4,4,1,1))
      plot(prioritisation$average_path, prioritisation$betweenness, pch=16,
           col = colorRampPalette(c("royalblue", "orange"))(nrow(full_site_list)-1)[1:(nrow(full_site_list)-1)],
           xlab="Average path length",
           ylab="Betweenness")

    }
  }

  return(list(prioritisation = prioritisation,
              network = full_network,
              AUC=AUC,
              site_list = full_site_list))

}
