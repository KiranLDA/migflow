#' Generate a random network
#'
#' @description This function generates a random migratory network
#'
#'
#' @param pop population size flowing through network
#' @param mean_dist average distance the species can travel (this can be estimated from real data)
#' @param sd_dist standard deviation of the distances a species can travel (this can be estimated from real data)
#' @param toplot TRUE/FALSE to determine whether the output is plotted or not
#' @param nbreeding Number of breeding sites
#' @param nwintering Number of wintering sites
#' @param nstop Number of stopover sites (shared during north and south migration)
#'
#' @return a list containting the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @import igraph
#' @importFrom stats rnorm runif
#' @export
randomDIRECTED <- function(nbreeding = 3,
                      nwintering = 4,
                      nstop = 10,
                      pop = 100000,
                      mean_dist = 1000,
                      sd_dist = 200,
                      toplot = TRUE){

  # # For testing purposes
  # nbreeding = 3
  # nwintering = 3
  # nstop = 10
  # pop = 100000
  # toplot= TRUE
  # # minforks = 2
  # # maxforks = 5
  # # mean_dist = 0
  # # sd_dist = 2000



  nsites = nbreeding+nwintering+nstop


  # generate random tracks
  tracks <- abs(rnorm(1000, mean_dist, sd_dist))


  # Create a fake list of sites where animals were seen at, with latitude, longitude and number of anumals seen there
  site_list <- data.frame(Lat= runif(nsites, min=-20, max=40),
                          Lon= runif(nsites, min=-5, max=5),
                          Pop=runif(nsites, min=500, max=10000),
                          B= 0,
                          SM=0,
                          NM=0,
                          NB=0)

  #sort according to latitude
  site_list = site_list[order(site_list$Lat, decreasing=T),]
  site_list$Site= 1:nsites

  site_list$B[1:nbreeding] = 1
  site_list$NB[(nrow(site_list)-nwintering+1):nrow(site_list)] = 1
  site_list$SM[site_list$B==0 & site_list$NB==0] = 1
  site_list$NM = site_list$SM

  # # add a dummy breeding and wintering site
  # site_list<- rbind(
  #   c(41,0,pop,0,0,0,0,0),
  #   site_list,
  #   c(-21,0,pop,0,0,0,0,9999))


  #-----------------------------
  #  South migration B -> NB
  #-----------------------------


  sites = site_list
  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)
  # Dist_P <- (max(dist)-dist) / max(dist)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )+0.01

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  network <- directedNET(network, include_diagonal = TRUE)

  # Ensure that nodes only flow into the next 1 or two neighbouring nodes
  for (i in 1:nrow(network)){
    idx = (i+1) : (i+1 + ceiling(runif(1,0,2))-1)
    idx = which(!(1:nrow(network) %in% idx))
    network[i, idx] = 0
  }

  SMnet <- t(apply(network,1,
                   function(x) x[which(!is.na(x))]/
                     sum(x,na.rm=TRUE)))



  #-----------------------------
  #  North migration NB -> B
  #-----------------------------

  sites = site_list[order(site_list$Lat, decreasing=FALSE),]

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )+0.01

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  network <- directedNET(network, include_diagonal = TRUE)

  # Ensure that nodes only flow into the next 1 or two neighbouring nodes
  for (i in 1:nrow(network)){
    idx = (i+1) : (i+1 + ceiling(runif(1,0,2))-1)
    idx = which(!(1:nrow(network) %in% idx))
    network[i, idx] = 0
  }

  NMnet <- t(apply(network,1,
                   function(x) x[which(!is.na(x))]/
                     sum(x,na.rm=TRUE)))



  #----------------------------------------
  # Add supersource and sink nodes
  #----------------------------------------

  #South
  SMnet <- addSUPERNODE(SMnet,
                        sources= site_list$Site[site_list$B ==1],
                        sinks = site_list$Site[site_list$NB ==1])

  index = as.numeric(names(which(SMnet["supersource",] == Inf)))
  SMnet["supersource",which(SMnet["supersource",] == Inf)] = site_list$Pop[index]/sum(site_list$Pop[index])
  index = as.numeric(names(which(SMnet[,"supersink"] == Inf)))
  SMnet[which(SMnet[,"supersink"] == Inf), "supersink"] = site_list$Pop[index]/sum(site_list$Pop[index])

  #North
  NMnet <- addSUPERNODE(NMnet,
                        sources= site_list$Site[site_list$NB ==1],
                        sinks = site_list$Site[site_list$B ==1])

  index = as.numeric(names(which(NMnet["supersource",] == Inf)))
  NMnet["supersource",which(NMnet["supersource",] == Inf)] = site_list$Pop[index]/sum(site_list$Pop[index])
  index = as.numeric(names(which(NMnet[,"supersink"] == Inf)))
  NMnet[which(NMnet[,"supersink"] == Inf), "supersink"] = site_list$Pop[index]/sum(site_list$Pop[index])



  colnames(NMnet)[1] = rownames(NMnet)[1] = "supersink"
  colnames(NMnet)[length(NMnet[1,])] = rownames(NMnet)[length(NMnet[1,])] = "supersource"


  site_list<- rbind(
    c(-21,0,pop,0,0,0,0,"supersink"),
    site_list,
    c(41,0,pop,0,0,0,0,"supersource"))

  # NMnet <- network

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

  network["Ssupersink","NB"]<- 1
  network["NB","Nsupersink"]<- 1




  network = network*pop

  #----------------------------------

  sites <- site_list


  weight <- graph_from_adjacency_matrix(network,  mode="directed", weighted = TRUE)

  # run the population through the network a forst time
  flow = max_flow(weight, source = V(weight)["Ssupersource"],
                  target = V(weight)["Nsupersource"], capacity = E(weight)$weight)

  neti <- weight
  E(neti)$weight <- flow$flow
  network <- as.matrix(as_adjacency_matrix(neti, attr="weight"))

  #created a weigted igraph network
  if (toplot == TRUE){
    # plot flow network
    nodes = get.edgelist(weight, names=TRUE)
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
    par(mfrow=c(1,1))
    par(mar=c(4,4,4,4))
    index=2:(nrow(sites)-1)
    plot(sites$Lon[index], sites$Lat[index], pch=16,
         cex=0)

    index=1:nrow(nodes)#2:(nrow(nodes)-1)#
    # Arrows(x0 = nodes$Lon_from[index],
    #        y0 = nodes$Lat_from[index],
    #        x1 = nodes$Lon_to[index],
    #        y1 = nodes$Lat_to[index],
    #        col= adjustcolor("royalblue3", alpha.f =  0.9))#,
    #        # lwd=(nodes$flow[index]/(max(nodes$flow)))*30)#,#/500
    segments(x0 = nodes$Lon_from[index],
             y0 = nodes$Lat_from[index],
             x1 = nodes$Lon_to[index],
             y1 = nodes$Lat_to[index],
             col= "black",#adjustcolor("royalblue3", alpha.f = 0.9),
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
    index=as.numeric(nodeflowplot$Category)+1
    points(sites$Lon[index],
           sites$Lat[index],
           pch=21,
           cex=(((nodeflowplot$x)/
                   as.numeric(max(nodeflowplot$x)))+0.4)*4,
           bg="orange", col="black")
  }

  return(list( network = network,
               tracks = tracks,
               sites = site_list  ))

}
