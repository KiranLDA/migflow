#' Prioritise sites according to population flow
#'
#' @description Prioritise sites according to population flow
#'
#' @param network adjacency network describing where birds are going
#' @param sites network sites including latitude and longitude, and population count at each site
#' @param returnNET TRUE/FALSE whether or not the network is a return network or not
#' @param toplot TRUE/FALSE to determine whether the output is plotted or not
#' @param pop number of individuals in the population
#'
#' @return a list containting the prioritisation (as a list of the sites in the order they are removed), the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#' pop=100000
#' rand_net = randomNET(nsites=15,pop=pop)
#'
#' # priotise sites according to flow through network
#' prioritiseFLOW(rand_net$network, rand_net$sites, returnNET=FALSE)
#'
#' # with a return network
#' net <- randomSTAR(nwintering=40, nbreeding=5,anadromous=FALSE)
#' network <- shortenNET(net$network, from = "Ssupersink", to = "Nsupersink")
#'
#'
#' colnames(network)[1]<- rownames(network)[1] <- "Ssupersource"
#' colnames(network)[length(colnames(network))] <-  "Nsupersink"
#' rownames(network)[length(rownames(network))] <- "Nsupersink"
#' sites <- net$sites
#'
#' # priotise sites according to flow through network
#' prioritiseFLOW(network, sites, toplot=TRUE, returnNET=TRUE)
#'
#'
#'
#' # with a return network
#' net <- randomPARALLEL(nwintering=40, sd_dist=2000, nbreeding=5)
#' network <- shortenNET(net$network, from = "Ssupersink", to = "Nsupersink")
#'
#'
#' colnames(network)[1]<- rownames(network)[1] <- "Ssupersource"
#' colnames(network)[length(colnames(network))] <-  "Nsupersink"
#' rownames(network)[length(rownames(network))] <- "Nsupersink"
#' sites <- net$sites
#'
#' # priotise sites according to flow through network
#' prioritiseFLOW(network, sites, toplot=TRUE, returnNET=TRUE)
#'
#'
#' # with a return network
#' net <- randomDIRECTED(nwintering=40, sd_dist=2000, nbreeding=5)
#' network <- shortenNET(net$network, from = "Ssupersink", to = "Nsupersink")
#'
#' colnames(network)[1]<- rownames(network)[1] <- "Ssupersource"
#' colnames(network)[length(colnames(network))] <-  "Nsupersink"
#' rownames(network)[length(rownames(network))] <- "Nsupersink"
#' sites <- net$sites
#'
#' # priotise sites according to flow through network
#' prioritiseFLOW(network, sites, toplot=TRUE, returnNET=TRUE)
#'
#' @import igraph
#' @importFrom graphics par segments points
#' @importFrom stats aggregate
#' @importFrom grDevices colorRampPalette
#'
#' @export
prioritiseFLOW <- function(network,
                           sites,
                           # method ="igraph",
                           toplot = TRUE,
                           pop = 100000,
                           # capacity,
                           returnNET = TRUE){
  # if (method == "igraph"){
  if(returnNET==FALSE){
    #created a weigted igraph network
    weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

    if(!exists("capacity")){
      capacity <- E(weight)$weight
    }

    # run the population through the network a forst time
    flow = max_flow(weight, source = V(weight)["supersource"],
                    target = V(weight)["supersink"],
                    capacity = capacity)
                      #E(weight)$weight )

    # plot flow network
    if (toplot == TRUE){
      par(mfrow=c(1,3))
      par(mar=c(0,0,0,0))
      index=2:(nrow(sites)-1)
      plot(sites$Lon[index], sites$Lat[index], pch=16,
           cex=0,
           axes=F)
    }

    nodes = as_edgelist(weight, names=TRUE)
    nodes = as.data.frame(nodes)
    nodes$flow = flow$flow

    nodes$Lat_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lat[sites$Site %in% nodes[i,1]])))
    nodes$Lon_from = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lon[sites$Site %in% nodes[i,1]])))
    nodes$Lat_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lat[sites$Site %in% nodes[i,2]])))
    nodes$Lon_to   = unlist(lapply(1:nrow(nodes), function(i) as.numeric(sites$Lon[sites$Site %in% nodes[i,2]])))

    nodes2 = nodes[nodes$V1 != "supersource" & nodes$V2 != "supersink" ,]

    if (toplot == TRUE){
      index=2:(nrow(nodes2)-1)
      segments(x0 = nodes2$Lon_from[index],
               y0 = nodes2$Lat_from[index],
               x1 = nodes2$Lon_to[index],
               y1 = nodes2$Lat_to[index],
               lwd=(nodes2$flow[index]/pop)*30)
    }

    # sort sites by flow
    nodeflow = merge(aggregate(nodes$flow, by=list(Category=as.character(nodes$V1)), FUN=sum),
                     aggregate(nodes$flow, by=list(Category=as.character(nodes$V2)), FUN=sum), all=T)
    nodeflow$x = as.numeric(nodeflow$x)
    nodeflow = data.frame( unique(as.matrix(nodeflow[ , 1:2 ]) ))
    nodeflow$x = as.numeric(as.character(nodeflow$x))
    nodeflow = nodeflow[nodeflow$Category != "supersource" & nodeflow$Category != "supersink",]

    # specify the site to remove, in this case, it is the site which contributes less to population flow
    to_remove = which(nodeflow$x %in% min(nodeflow$x))

    # empty dataset to store output
    prioritisation <- data.frame(Site=as.character(nodeflow$Category[to_remove]),
                                 Pop_Flow =   flow$value,
                                 Site_Flow =   nodeflow$x[to_remove])
    # make sure it is numeric
    nodeflow$Category = as.numeric(as.character(nodeflow$Category))

    # plot sites
    if (toplot == TRUE){
      nodeflowplot = nodeflow[order(nodeflow$Category),]
      index=as.numeric(nodeflowplot$Category)+1
      points(sites$Lon[index],
             sites$Lat[index],
             pch=21,
             cex=(((nodeflowplot$x)/
                     as.numeric(max(nodeflowplot$x)))+0.4)*4,
             bg="orange", col="black")
    }

    # store full network for results, before things get removed in the optimisation
    full_network = network
    full_site_list = sites

    # prioritisation
    net_remove = which(colnames(network) %in% nodeflow$Category[to_remove])
    network = network[-net_remove,-net_remove ]
    sites = sites[!(sites$Site %in% as.character(nodeflow$Category[to_remove])),]


    while(nrow(sites) > 2 & sum(network)>0){


      weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

      flow = max_flow(weight, source = V(weight)["supersource"],
                      target = V(weight)["supersink"], capacity = E(weight)$weight )

      nodes = as_edgelist(weight, names=TRUE)
      nodes = as.data.frame(nodes)
      nodes$flow = flow$flow

      nodeflow = merge(aggregate(nodes$flow, by=list(Category=as.character(nodes$V1)), FUN=sum),
                       aggregate(nodes$flow, by=list(Category=as.character(nodes$V2)), FUN=sum), all=T)
      nodeflow$x = as.numeric(nodeflow$x)

      nodeflow = data.frame( unique(as.matrix(nodeflow[ , 1:2 ]) ))
      nodeflow$x = as.numeric(as.character(nodeflow$x))

      nodeflow = nodeflow[nodeflow$Category != "supersource" & nodeflow$Category != "supersink",]
      to_remove = which(nodeflow$x %in% min(nodeflow$x))

      prioritisation <- rbind(prioritisation,
                              data.frame(Site=as.character(nodeflow$Category[to_remove]),
                                         Pop_Flow =   flow$value,
                                         Site_Flow =   nodeflow$x[to_remove]))

      net_remove = which(colnames(network) %in% nodeflow$Category[to_remove])
      network[-net_remove,-net_remove ]
      network = network[-net_remove,-net_remove ]
      sites = sites[!(sites$Site %in% as.character(nodeflow$Category[to_remove])),]

    }
    y = prioritisation$Pop_Flow/prioritisation$Pop_Flow[1]
    x = (1:length(prioritisation$Pop_Flow))/length(prioritisation$Pop_Flow)
    AUC = curvCALC(x,y)

    if (toplot == TRUE){
      par(mar=c(4,4,1,1))
      plot(prioritisation$Pop_Flow, type="o",pch=16,
           xlab="# Sites Removed",
           ylab="Population size",
           main = AUC)

      par(mar=c(4,4,1,1))
      plot(prioritisation$Site_Flow, c(0,abs(diff(prioritisation$Pop_Flow))), pch=16,
           col = colorRampPalette(c("royalblue", "orange"))(nrow(full_site_list)-1)[1:(nrow(full_site_list)-1)],
           xlab="Site carrying capacity",
           ylab="Change in population size when site lost")
    }
  }


  if(returnNET==TRUE){
    #created a weigted igraph network
    # network[network==0]<-0.01
    weight <- graph_from_adjacency_matrix(network, mode="directed", weighted = TRUE)


#     if(!exists("capacity")){
#       capacity <- E(weight)$weight
#     }

    sourcename = rownames(network)[1]
    sinkname =  colnames(network)[length(colnames(network))]

    # run the population through the network a forst time
    flow = max_flow(weight, source = V(weight)[sourcename],
                    target = V(weight)[sinkname],
                    capacity = E(weight)$weight )
                   # capacity = capacity)




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
    to_remove = which(nodeflow$x %in% min(nodeflow$x))

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




    # prioritisation

    net_remove = which(colnames(network) %in% as.character(nodeflow$Category[to_remove]))
    network = network[-net_remove,-net_remove ]
    sites = sites[!(sites$Site %in% as.character(nodeflow$Category[to_remove])),]


    while(nrow(sites) > 2 & sum(network)>0){

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
      to_remove = which(nodeflow$x %in% min(nodeflow$x))

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

      net_remove = which(colnames(network) %in% as.character(nodeflow$Category[to_remove]))
      # network[-net_remove,-net_remove ]
      network = network[-net_remove,-net_remove ]
      sites = sites[!(sites$Site %in% as.character(nodeflow$Category[to_remove])),]

    }

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
  return(list(# method = method,
               prioritisation = prioritisation,
               network = full_network,
               AUC=AUC,
               site_list = full_site_list  ))
}
