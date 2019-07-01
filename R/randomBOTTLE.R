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
randomBOTTLE <- function(nbreeding = 10,
                           nwintering = 10,
                           nstop = 30,
                           pop = 100000,
                           mean_dist = 0,
                           sd_dist = 500,
                           toplot = TRUE){

  # # # For testing purposes
  # nbreeding = 10
  # nwintering = 10
  # nstop = 30
  # pop = 100000
  # mean_dist = 0
  # sd_dist = 800
  # toplot = TRUE


  nsites = nbreeding+nwintering+nstop


  # generate random tracks
  tracks <- abs(rnorm(1000, mean_dist, sd_dist)) # rgamma(1000,2.3,0.005)#
  # hist(tracks)

  B <- data.frame(Lat = runif(nbreeding, min=30, max=40),
                  Lon= runif(nbreeding, min=-20, max=20),
                  Pop=runif(nbreeding, min=500, max=10000),
                  B = 1,
                  SM = 0,
                  NB = 0,
                  NM = 0)


  NB <- data.frame(Lat = runif(nwintering, min=-20, max=-10),
                   Lon= runif(nwintering, min=-20, max=20),
                   Pop=runif(nwintering, min=500, max=10000),
                   B = 0,
                   SM = 0,
                   NB = 1,
                   NM = 0)

  STP <- data.frame(Lat = runif(nstop, min=-10, max=30),
                    Lon= runif(nstop, min=-20, max=20),
                    Pop=runif(nstop, min=500, max=10000),
                    B = 0,
                    SM = 1,
                    NB = 0,
                    NM = 1)

  sites = rbind(B, STP, NB)
  sites = sites[order(sites$Lat, decreasing=T),]
  sites$Site= 1:nsites
  site_list = sites


  #Add a bottleneck site
  bottleneck <- round(runif(1, 1, nstop))
  site_list$Pop[site_list$SM==1][bottleneck] <- pop


  # sort according to latitude
  site_list = site_list[order(site_list$Lat, decreasing=T),]
  site_list$Site= 1:nsites

  site_list$B[1:nbreeding] = 1
  site_list$NB[(nrow(site_list)-nwintering+1):nrow(site_list)] = 1
  site_list$SM[site_list$B==0 & site_list$NB==0] = 1
  site_list$NM = site_list$SM


  #-----------------------------
  #  South migration B -> NB
  #-----------------------------


  sites = site_list

  # different!
  bottleneck_site <- site_list$Site[site_list$SM==1][bottleneck]
  bottleneck_idx <- which(sites$Site == bottleneck_site)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F) +0.01
  # different!
  Dist_P[,bottleneck_idx] <- (max(dist)-dist[,bottleneck_idx]) / max(dist)
  Dist_P[bottleneck_idx,] <- (max(dist)-dist[bottleneck_idx,]) / max(dist)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)
  # different!
  Pop_P[,bottleneck_idx] <- 1
  Pop_P[bottleneck_idx,] <- 1

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )#+0.01
  # different!
  # Azi_P[,bottleneck_idx] <- 1
  # Azi_P[bottleneck_idx,] <- 1

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <-  Azi_P  * Pop_P * Dist_P

  # different!
  # Dist_P <- (max(dist)-dist) / max(dist)
  # network[, which(sites$SM == 1 & sites$Pop == pop)] <- Azi_P[, which(sites$SM == 1 & sites$Pop == pop)]*Pop_P[, which(sites$SM == 1 & sites$Pop == pop)]*Dist_P[, which(sites$SM == 1 & sites$Pop == pop)]
  # network[which(sites$SM == 1 & sites$Pop == pop), ] <- Azi_P[which(sites$SM == 1 & sites$Pop == pop), ]*Pop_P[which(sites$SM == 1 & sites$Pop == pop), ]*Dist_P[which(sites$SM == 1 & sites$Pop == pop), ]


  # Make the network directed
  network <- directedNET(network, include_diagonal = TRUE)

  SMnet <- t(apply(network,1,
                   function(x) x[which(!is.na(x))]/
                     sum(x,na.rm=TRUE)))



  #-----------------------------
  #  North migration NB -> B
  #-----------------------------

  sites = site_list[order(site_list$Lat, decreasing=FALSE),]

  # different!
  bottleneck_idx <- which(sites$Site == bottleneck_site)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F) +0.01
  # different!
  Dist_P[,bottleneck_idx] <- (max(dist)-dist[,bottleneck_idx]) / max(dist)
  Dist_P[bottleneck_idx,] <- (max(dist)-dist[bottleneck_idx,]) / max(dist)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)
  # different!
  Pop_P[,bottleneck_idx] <- 1
  Pop_P[bottleneck_idx,] <- 1

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )#+0.01
  # # different!
  # Azi_P[,bottleneck_idx] <- 1
  # Azi_P[bottleneck_idx,] <- 1

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <-  Azi_P  *Pop_P * Dist_P

  # Make the network directed
  network <- directedNET(network, include_diagonal = TRUE)

  # different!
  # Dist_P <- (max(dist)-dist) / max(dist)
  # network[, which(sites$SM == 1 & sites$Pop == pop)] <- Azi_P[, which(sites$SM == 1 & sites$Pop == pop)]*Pop_P[, which(sites$SM == 1 & sites$Pop == pop)] * Dist_P[, which(sites$SM == 1 & sites$Pop == pop)]
  # network[which(sites$SM == 1 & sites$Pop == pop), ] <- Azi_P[which(sites$SM == 1 & sites$Pop == pop), ]*Pop_P[which(sites$SM == 1 & sites$Pop == pop), ]*Dist_P[which(sites$SM == 1 & sites$Pop == pop), ]

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
  SMnet["supersource",which(SMnet["supersource",] == Inf)] = site_list$Pop[index]/sum(site_list$Pop[index]) #pop / nbreeding#
  index = as.numeric(names(which(SMnet[,"supersink"] == Inf)))
  SMnet[which(SMnet[,"supersink"] == Inf), "supersink"] = site_list$Pop[index]/sum(site_list$Pop[index]) # pop / nwintering#

  #North
  NMnet <- addSUPERNODE(NMnet,
                        sources= site_list$Site[site_list$NB ==1],
                        sinks = site_list$Site[site_list$B ==1])

  index = as.numeric(names(which(NMnet["supersource",] == Inf)))
  NMnet["supersource",which(NMnet["supersource",] == Inf)] = site_list$Pop[index]/sum(site_list$Pop[index]) #pop / nwintering#
  index = as.numeric(names(which(NMnet[,"supersink"] == Inf)))
  NMnet[which(NMnet[,"supersink"] == Inf), "supersink"] = site_list$Pop[index]/sum(site_list$Pop[index]) #pop / nbreeding#



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

  # neti <- weight
  # E(neti)$weight <- flow$flow
  # network <- as.matrix(as_adjacency_matrix(neti, attr="weight"))

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
         cex=0, xlab="", ylab="", xaxt="n", yaxt = "n",
         frame.plot=FALSE)

    index=1:nrow(nodes)#2:(nrow(nodes)-1)#
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
    colorz = ifelse(sites$B[index]==1,"royalblue",ifelse(sites$NB[index]==1,"orange","gray"))
    points(sites$Lon[index],
           sites$Lat[index],
           pch=21,
           cex=(((nodeflowplot$x)/
                   as.numeric(max(nodeflowplot$x)))+0.4)*4,
           bg=colorz , col="black")
  }

  return(list( network = network,
               tracks = tracks,
               sites = site_list  ))

}
