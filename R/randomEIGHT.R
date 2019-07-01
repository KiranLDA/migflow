#' Generate a random network
#'
#' @description This function generates a random migratory network
#'
#' @param nbreeding number of breeding sites.
#' @param nwintering number of nonbreeding residency/wintering sites.
#' @param npre number of sites used duting pre-breeding (north) migration.
#' @param npost number of sites used duting post-breeding (south) migration.
#' @param toplot TRUE/FALSE to determine whether the output is plotted or not
#' @param pop population size flowing through network
#' @param mean_dist mean distance an animal can travel
#' @param sd_dist standard deviation of the distance an animal can travel
#'
#'
#' @return a list containting the network which was randomly generated,
#' the tracks that were randomly generated, and the sites that were randomly generated for animals to use.
#'
#' @examples
#' par(mfrow=c(1,1))
#' randomLOOP(toplot=TRUE)
#' randomLOOP(nbreeding = 3,
#'            nwintering = 3,
#'            npre = 3,
#'            npost = 3,
#'            toplot=TRUE)
#'
#'
#' @import igraph
#' @importFrom stats runif
#' @export
randomEIGHT <- function(pop = 100000,
                        toplot= TRUE,
                        nbreeding = 5,
                        nwintering = 10,
                        npre = 20,
                        npost = 20,
                        mean_dist = 2000,
                        sd_dist =1000){

  # For testing purposes
  # pop = 100000
  # toplot= TRUE
  # nbreeding = 8
  # nwintering = 10
  # npre = 18
  # npost = 18
  # minforks = 2
  # maxforks = 5
  #
  # mean_dist = 0
  # sd_dist = 5000

  # generate random tracks
  tracks <- abs(rnorm(100, mean_dist, sd_dist))


  B <- data.frame(Lat = runif(nbreeding, min=30, max=40),
                  Lon= runif(nbreeding, min=-20, max=20),
                  Pop=runif(nbreeding, min=500, max=10000),
                  B = 1,
                  SM = 0,
                  NB = 0,
                  NM = 0)
  # B <- B[order(B$Lon, decreasing=FALSE),]
  B <- B[order(B$Lat, decreasing=TRUE),]
  B$Site = 1: nbreeding

  SM1 <- data.frame(Lat = runif(ceiling(npost/2), min=10, max=30),
                    Lon = runif(ceiling(npost/2), min=10, max=20),
                    Pop = runif(ceiling(npost/2), min=500, max=10000),
                    B = 0,
                    SM = 1,
                    NB = 0,
                    NM = 0)
  SM1 <- SM1[order(SM1$Lat, decreasing=TRUE),]
  SM1$Site = (nbreeding + 1): (nbreeding + ceiling(npost/2))

  SM2 <- data.frame(Lat = runif(floor(npost/2), min=-10, max=10),
                    Lon = runif(floor(npost/2), min=-20, max=-10),
                    Pop = runif(floor(npost/2), min=500, max=10000),
                    B = 0,
                    SM = 1,
                    NB = 0,
                    NM = 0)
  SM2 <- SM2[order(SM2$Lat, decreasing=TRUE),]
  SM2$Site = (nbreeding + ceiling(npost/2) + 1): (nbreeding + ceiling(npost/2) + floor(npost/2) )

  NM1 <- data.frame(Lat = runif(floor(npre/2), min=-10, max=10),
                    Lon = runif(floor(npre/2), min=10, max=20),
                    Pop = runif(floor(npre/2), min=500, max=10000),
                    B = 0,
                    SM = 0,
                    NB = 0,
                    NM = 1)

  NM1 <- NM1[order(NM1$Lat, decreasing=FALSE),]
  NM1$Site = (nbreeding + ceiling(npost/2) + floor(npost/2) + 1) : (nbreeding + ceiling(npost/2) +
                                                                      floor(npost/2) + floor(npre/2) )


  NM2 <- data.frame(Lat = runif(ceiling(npre/2), min=10, max=30),
                    Lon = runif(ceiling(npre/2), min=-20, max=-10),
                    Pop = runif(ceiling(npre/2), min=500, max=10000),
                    B = 0,
                    SM = 0,
                    NB = 0,
                    NM = 1)

  NM2 <- NM2[order(NM2$Lat, decreasing=FALSE),]
  NM2$Site = (nbreeding + ceiling(npost/2) +
                floor(npost/2) + floor(npre/2) + 1) : (nbreeding + ceiling(npost/2) + floor(npost/2) +
                                                         floor(npre/2) + ceiling(npre/2))

  NB <- data.frame(Lat = runif(nwintering, min=-20, max=-10),
                   Lon= runif(nwintering, min=-20, max=20),
                   Pop=runif(nwintering, min=500, max=10000),
                   B = 0,
                   SM = 0,
                   NB = 1,
                   NM = 0)
  # NB <- NB[order(NB$Lon, decreasing=FALSE),]
  NB <- NB[order(NB$Lat, decreasing=FALSE),]
  NB$Site = (nbreeding + ceiling(npost/2) +
               floor(npost/2) + floor(npre/2) +
               ceiling(npre/2) + 1 ) : (nbreeding + ceiling(npost/2) +
                                          floor(npost/2) + floor(npre/2) +
                                          ceiling(npre/2) + nwintering )


  sites = rbind(B, SM1, SM2, NM1, NM2, NB)
  site_list = sites

  #--------------------------------------------------------
  #            SOUTH MIGRATION
  #--------------------------------------------------------

  #----------------------------------
  # Step 1
  #----------------------------------

  sites <- rbind( B, SM1)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  network1 <- directedNET(network, include_diagonal = TRUE)
  network1 <- cbind(network1,
                    matrix(NA, nrow(network1), length(SM2$Site)),
                    matrix(NA, nrow(network1), length(NB$Site)))
  network1 <- rbind(network1,
                    matrix(NA, length(SM2$Site), ncol(network1)),
                    matrix(NA, length(NB$Site), ncol(network1)))
  colnames(network1) = rownames(network1) = c(B$Site,SM1$Site, SM2$Site, NB$Site)
  # for(i in SM2$Site) { network1[,i] <- 0}
  # for(i in SM2$Site) { network1[i,] <- 0}

  #----------------------------------
  # Step 2
  #----------------------------------

  sites <- rbind( SM1, SM2)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  # network2 <- directedNET(network, include_diagonal = TRUE)
  network2 <- directedNET(network, include_diagonal = TRUE)
  network2 <- cbind(matrix(NA, nrow(network2), length(B$Site)),
                    network2,
                    matrix(NA, nrow(network2), length(NB$Site))
  )
  network2 <- rbind(matrix(NA, length(B$Site) , ncol(network2)),
                    network2,
                    matrix(NA, length(NB$Site), ncol(network2)))
  colnames(network2) = rownames(network2) = c(B$Site,SM1$Site, SM2$Site, NB$Site)



  #----------------------------------
  # Step 3
  #----------------------------------

  sites <- rbind( SM2, NB)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  # network2 <- directedNET(network, include_diagonal = TRUE)
  network3 <- directedNET(network, include_diagonal = TRUE)
  network3 <- cbind(matrix(NA, nrow(network3), length(B$Site)),
                    matrix(NA, nrow(network3), length(SM1$Site)),
                    network3)
  network3 <- rbind(matrix(NA, length(B$Site) , ncol(network3)),
                    matrix(NA, length(SM1$Site), ncol(network3)),
                    network3)
  colnames(network3) = rownames(network3) = c(B$Site,SM1$Site, SM2$Site, NB$Site)


  SMnet <- pmax(network1, network2, na.rm = TRUE)
  SMnet <- pmax(SMnet, network3, na.rm = TRUE)
  SMnet[is.na(SMnet)] <-0

  SMnet <- t(apply(SMnet,1,
                   function(x) x[which(!is.na(x))]/
                     sum(x,na.rm=TRUE)))


  #----------------------------------------
  # Add supersource and sink nodes
  #----------------------------------------

  SMnet <- addSUPERNODE(SMnet,
                        sources= site_list$Site[site_list$B ==1],
                        sinks = site_list$Site[site_list$NB ==1])

  index = as.numeric(names(which(SMnet["supersource",] == Inf)))
  SMnet["supersource",which(SMnet["supersource",] == Inf)] = site_list$Pop[index]/sum(site_list$Pop[index])
  index = as.numeric(names(which(SMnet[,"supersink"] == Inf)))
  SMnet[which(SMnet[,"supersink"] == Inf), "supersink"] = site_list$Pop[index]/sum(site_list$Pop[index])

  #--------------------------------------------------------
  #            NORTH MIGRATION
  #--------------------------------------------------------

  #----------------------------------
  # Step 1
  #----------------------------------

  sites <- rbind( NB, NM1)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  network1 <- directedNET(network, include_diagonal = TRUE)
  network1 <- cbind(network1,
                    matrix(NA, nrow(network1), length(NM2$Site)),
                    matrix(NA, nrow(network1), length(B$Site)))
  network1 <- rbind(network1,
                    matrix(NA, length(NM2$Site), ncol(network1)),
                    matrix(NA, length(B$Site), ncol(network1)))
  colnames(network1) = rownames(network1) = c(NB$Site,NM1$Site, NM2$Site, B$Site)
  # for(i in SM2$Site) { network1[,i] <- 0}
  # for(i in SM2$Site) { network1[i,] <- 0}

  #----------------------------------
  # Step 2
  #----------------------------------

  sites <- rbind( NM1, NM2)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P * Pop_P * Azi_P

  # Make the network directed
  # network2 <- directedNET(network, include_diagonal = TRUE)
  network2 <- directedNET(network, include_diagonal = TRUE)
  network2 <- cbind(matrix(NA, nrow(network2), length(NB$Site)),
                    network2,
                    matrix(NA, nrow(network2), length(B$Site))
  )
  network2 <- rbind(matrix(NA, length(NB$Site) , ncol(network2)),
                    network2,
                    matrix(NA, length(B$Site), ncol(network2)))
  colnames(network2) = rownames(network2) = c(NB$Site,NM1$Site, NM2$Site, B$Site)



  #----------------------------------
  # Step 3
  #----------------------------------

  sites <- rbind( NM2, B)

  # create a distance matrix based on these data
  dist <- point2DIST(sites)

  # calculate the probability of going between these sites given the distance the animal can travel
  Dist_P <- distPROB(tracks, dist, adjust=1, plot=F)

  # Calculate prioritisation of population using a site
  Pop_P <- nodePopPROP(sites, population = pop)

  #Calculate the azimuth angle
  Azi_P <- absAZIMUTH(dist, lonlats=sites )

  # make birds/animals prefer sites which a larger prioritisation of the population has been seen and where the distance is better
  network <- Dist_P  * Azi_P #* Pop_P


  # Make the network directed
  # network2 <- directedNET(network, include_diagonal = TRUE)
  network3 <- directedNET(network, include_diagonal = TRUE)

  network3 <- cbind(matrix(NA, nrow(network3), length(NB$Site)),
                    matrix(NA, nrow(network3), length(NM1$Site)),
                    network3)
  network3 <- rbind(matrix(NA, length(NB$Site) , ncol(network3)),
                    matrix(NA, length(NM1$Site), ncol(network3)),
                    network3)
  colnames(network3) = rownames(network3) = c(NB$Site,NM1$Site, NM2$Site, B$Site)


  NMnet <- pmax(network1, network2, na.rm = TRUE)
  NMnet <- pmax(NMnet, network3, na.rm = TRUE)
  NMnet[is.na(NMnet)] <-0
  NMnet <- t(apply(NMnet,1,
                   function(x) x[which(!is.na(x))]/
                     sum(x,na.rm=TRUE)))

  #----------------------------------------
  # Add supersource and sink nodes
  #----------------------------------------


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
#
#   neti <- weight
#   E(neti)$weight <- flow$flow
#   network <- as.matrix(as_adjacency_matrix(neti, attr="weight"))

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

