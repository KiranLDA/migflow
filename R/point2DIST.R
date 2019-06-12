#' Make a distance matrix
#'
#' @param dta This is a dataframe. Columns of interest must be named: Site, Lat (for latitude), Lon (for longitude).
#' @param miles TRUE/FALSE. whether or not to calculate distance in miles. if FALSE, then estimates in kilometres
#' @param proj4string sp projection for latitude longitude
#'
#' @return A distance matrix in kilometres between all pairs sites
#'
#' @examples
#' dta <- data.frame(Site= LETTERS[1:4], Lat= 1:4, Lon= 5:8, Pop=100:103)
#' point2DIST(dta)
#'
#' @importFrom fields rdist.earth
#' @importFrom sp SpatialPoints CRS
#'
#' @export
point2DIST <- function(dta, miles=FALSE,
                       proj4string = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "){
  # Project latitude and longitude
  pts  = SpatialPoints(coords = cbind(dta$Lat, dta$Lon),
                       proj4string=CRS(proj4string))
  # Calculate earth distances and spit out distance matrix
  dist = rdist.earth(pts@coords, miles=miles)
  colnames(dist) = rownames(dist) = dta$Site
  dist[is.na(dist)] = 0
  return(dist)
}
