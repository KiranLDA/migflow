#' Prioritise sites according to population flow
#'
#' @description Prioritise sites according to population flow
#'
#' @param network adjacency network describing where birds are going
#' @param sites network sites including latitude and longitude, and population count at each site
#' @param method max_flow method, currently that of igraph, but will also support gurobi or lpsolve
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
#' prioritiseFLOW(rand_net$network, rand_net$sites)
#'
#' # with fewer edges
#' pop=100000
#' rand_net = randomNET(nsites=15,pop=pop, nedges=40)
#'
#' # priotise sites according to flow through network
#' prioritiseFLOW(rand_net$network, rand_net$sites)
#'
#' @import igraph
#' @importFrom graphics par segments points
#' @importFrom stats aggregate
#'
#' @export
prioritiseFLOW <- function(network, sites, method ="igraph", toplot = TRUE, pop = 100000){
  if (method == "igraph"){

    #created a weigted igraph network
    weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

    # run the population through the network a forst time
    flow = max_flow(weight, source = V(weight)["supersource"],
                    target = V(weight)["supersink"], capacity = E(weight)$weight )

    # plot flow network
    if (toplot == TRUE){
      par(mfrow=c(1,2))
      par(mar=c(0,0,0,0))
      index=2:(nrow(sites)-1)
      plot(sites$Lon[index], sites$Lat[index], pch=16,
           cex=0,
           axes=F)
    }

    nodes = get.edgelist(weight, names=TRUE)
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

      nodes = get.edgelist(weight, names=TRUE)
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
    if (toplot == TRUE){
      par(mar=c(4,4,1,1))
      plot(prioritisation$Pop_Flow, type="o",pch=16,
           xlab="# Sites Removed",
           ylab="Population size")
    }
  }
  return(list( method = method,
               prioritisation = prioritisation,
               network = full_network,
               site_list = full_site_list  ))
}
