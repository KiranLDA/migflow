#' Calculate northness/southness of a site relative to another (from latitude longitude)
#'
#' @description Calculate the absolute azimuth between two sites
#'
#' @param dist This is a distance matrix between two sites based on lonlats
#' @param lonlats a dataframe containing columns entited Site, Lon, Lat, Pop
#'
#' @return A distance dist in kilometres between all pairs sites
#'
#' @examples
#' dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
#' dist <- point2DIST(dta)
#' absAZIMUTH(diat=dist, lonlats=dta)
#'
#' @importFrom maptools gzAzimuth
#'
#' @export
absAZIMUTH <- function(dist, lonlats){
 output=dist

 index0 <- cbind(rep(1:nrow(dist), ncol(dist)), rep(1:ncol(dist), each = nrow(dist)),
                rep(rownames(dist), ncol(dist)), rep(colnames(dist), each = nrow(dist)))

 index1 <- which(apply(index0, 1, function(x) !any(x%in%c("supersource","supersink"))))
 index0 <- cbind(index0,c(dist),c(dist))
 index0[index1,6] <- unlist(lapply(index1, function(x){
   xLon = lonlats$Lon[as.character(lonlats$Site) == index0[x,3]]
   xLat = lonlats$Lat[as.character(lonlats$Site) == index0[x,3]]
   yLon = lonlats$Lon[as.character(lonlats$Site) == index0[x,4]]
   yLat = lonlats$Lat[as.character(lonlats$Site) == index0[x,4]]
   z = gzAzimuth(c(xLon,xLat),c(yLon,yLat))
   if(is.na(z)) z=90
   return(abs(cos(z*pi/180)))}))

 output[which(dist==0)] = 0

 output[is.na(output)]= 0
 return(output)
}
