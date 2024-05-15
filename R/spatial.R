pullAndAverage <- function(ensNum,nulls,distWeights,areaWeights){
  #get null vals
  av <- map_dbl(nulls,ensNum)

  #calculate weights
  w <- distWeights + areaWeights

  w <- w/sum(w,na.rm = TRUE)

  return(sum(w * av,na.rm = TRUE))
}

gn2 <- function(x,type = "Above"){
  if(is.null(x[[paste0("pval",type)]])){
    return(NA)
  }else{
    if(is.na(x$pvalAbove)){
      return(NA)
    }else{
      return(getNulls(x$resultsUsed[[paste0("nullDetectionWithUncertainty",type)]],
                      nrows = nrow(x$resultsUsed),
                      n.ens = 1000,
                      weights = x$resultsUsed$weight))
    }
  }
}

gaspari_cohn <- function(r, rmax) {
  # This function takes two arguments:
  # r: a vector of distances from the point of interest to other points
  # rmax: the maximum distance over which the function is defined

  r <- abs(r)

  # Define constants
  a = 2/3 * rmax
  b = 4/3 * rmax

  # Calculate the function
  out = rep(0, length(r))
  ind1 = which(r <= a)
  out[ind1] = 1 - 7/6 * (r[ind1]/a)^2 + 3/6 * (r[ind1]/a)^3

  ind2 = which(r > a & r <= b)
  out[ind2] = 1/3 * (2 - (r[ind2]/a))^3

  ind3 = which(r > b)
  out[ind3] = 0

  return(out)
}


#' Distance from grid cells to TS objects
#'
#' @param wts a lipd TS object
#' @param grid.resolution grid cell size resolution, eg. 1 deg x 1 deg
#' @param radius.group radius of filter for each grid cell in km
#'
#' @return distances
#' @export
#'
distanceGrid <- function(wts=NULL,
                         lon.range = c(-179.5, 179.5),
                         lat.range = c(-89.5, 89.5),
                         grid.resolution = 1){

  #prep run
  allLon <- seq(min(lon.range),max(lon.range), grid.resolution)
  allLat <- seq(min(lat.range),max(lat.range), grid.resolution)
  grid <- expand.grid(allLon,allLat) |> setNames(c("lon","lat"))

  getAllDistances <- function(lat,lon,wts){
    list(lat = lat,
         lon = lon,
         dist = data.frame(TSid = wts$paleoData_TSid,
                      distance = t(geosphere::distm(c(lon, lat), cbind(wts$geo_longitude, wts$geo_latitude),fun = geosphere::distHaversine)/1000))
         )

  }


  allDistances <- purrr::map2(grid$lat,grid$lon,getAllDistances,wts,.progress = TRUE)

  return(allDistances)
}


#' Data Coverage Plot
#'
#' @param distance.grid as produced by distanceGrid()
#' @param color.breaks vector of breaks for coloring plot
#'
#' @return plot data
#' @export
#'
plotDataDensity <- function(distance.grid=NULL,
                            color.breaks = c(1,2,4,8,16)){

  plotData <- data.frame("TScount"=unlist(lapply(distance.grid, function(x) sum(!is.na(unique(x$dist$TSid))))),
                         "lat"=unlist(lapply(distance.grid, function(x) x$lat)),
                         "lon"=unlist(lapply(distance.grid, function(x) x$lon)))


  world <- map_data("world")


  p1 <- ggplot(data = plotData[plotData$TScount>0,],
               mapping = aes(x=lon,
                             y=lat,
                             fill=TScount)) +
    scale_fill_viridis_b(trans="log", breaks=color.breaks)+
    geom_tile() +
    geom_polygon(data=world, aes(x = long, y = lat, group = group), color="black", fill=alpha("white",0.1), inherit.aes = FALSE) +
    geom_point(inherit.aes = FALSE, data = wts, mapping = aes(y=geo_latitude,x=geo_longitude),color="white")+
    xlab("") +
    ylab("") +
    #title(main = "Data Coverage")+
    theme(legend.title = element_blank(),
          title = element_text("Data Coverage"))+
    coord_map(xlim=c(-180,180), ylim = c(-90,90), projection="mollweide")

  print(p1)

  returns <- list(plot = p1, plot.data=plotData)

  return(returns)
}




#' Title
#'
#' @param allNulls all the null data
#' @param nrows number of rows
#' @param n.ens number of ensembles
#' @param weights a vector to optionally weight the nulls
#'
#' @return a vector of null outputs
#' @export
getNulls <- function(allNulls,nrows,n.ens,weights = NA){

  if(all(is.na(weights))){
    weights <- rep(1,times = nrows)
  }

  if(length(allNulls) == 0){
stop("the null is missing")
  }

  weightMat <- rep(weights,times = length(allNulls[[1]])) %>%
    matrix(nrow = nrows, byrow = FALSE)

  nullMat <- matrix(unlist(allNulls),nrow = nrows,byrow = TRUE)

  nulls <- colSums(nullMat * weightMat)

  while(length(nulls) < n.ens){
    #shuffle column orders to increase sample size.
    shuffMat <- nullMat
    for(i in 1:nrow(nullMat)){
      shuffMat[i,] <- nullMat[i,sample.int(ncol(nullMat),replace = FALSE)]
    }
    nulls <- c(nulls,colSums(shuffMat * weightMat))
  }

  return(nulls)
}

#' Calculate the significance of multiple tests
#'
#' @param events a data.frame of events to consider, output of detectExcursion or detectShift
#' @param weights optionally weight the events
#' @param n.ens number of ensembles to consider
#'
#' @return data frame of results
#' @export
calculateMultiTestSignificance <- function(events,weights = NA,n.ens = 1000){

  #only include events with results
  events$weights <- weights
  events <- dplyr::filter(events,!is.na(event_probability))

  weights <- events$weights

  if(nrow(events) == 0){
    out <- list()
    out$pvalAbove <- NA
    out$pvalBelow <- NA
    out$pvalEither <- NA
    out$pvalBoth <- NA
    out$pvalNet <- NA
    return(out)
  }


  if(all(is.na(weights))){
    weights <- rep(1,times = nrow(events))
  }

  #make weights sum to 1
  weights <- weights/sum(weights)

  #get the sum of the event testing:
  allEventPos <- sum(events$event_probability_positive * weights,na.rm = FALSE)
  allEventNeg <- sum(events$event_probability_negative * weights,na.rm = FALSE)
  allEventEither <- sum(events$event_probability_either * weights,na.rm = FALSE)
  allEventBoth <- sum(events$event_probability_both * weights,na.rm = FALSE)
  allEventNet <- allEventPos - allEventNeg


  #compare the sum of the nulls
  #add in weights
  if(length(events$null_probability) == 0){
    stop("empty null")
  }

  nullsPos <- getNulls(events$null_probability_positive,nrow(events),n.ens,weights = weights)
  nullsNeg <- getNulls(events$null_probability_negative,nrow(events),n.ens,weights = weights)
  nullsEither <- getNulls(events$null_probability_either,nrow(events),n.ens,weights = weights)
  nullsBoth <- getNulls(events$null_probability_both,nrow(events),n.ens,weights = weights)
  nullsNet <- nullsPos - nullsNeg

  out <- list()
  out$pvalPos <- kdePval(nullsPos, allEventPos)$pval
  out$pvalNeg <- kdePval(nullsNeg, allEventNeg)$pval
  out$pvalEither <- kdePval(nullsEither, allEventEither)$pval
  out$pvalBoth <- kdePval(nullsBoth, allEventBoth)$pval
  out$pvalNet <- kdePval(nullsNet,allEventNet)$pval

  out$allEventPos <- allEventPos
  out$allEventNeg <- allEventNeg
  out$allEventBoth <- allEventBoth
  out$allEventEither <- allEventEither
  out$allEventNet <- allEventNet

  out$clPos <- quantile(nullsPos,probs = c(.975))
  out$clNeg <- quantile(nullsPos,probs = c(.975))
  out$clEither <- quantile(nullsEither,probs = c(.9,.95,.99))
  out$clBoth <- quantile(nullsBoth,probs = c(.9,.95,.99))
  out$clNet <- quantile(nullsNet,probs = c(.005,.025,.05,.95,.975,.995))


  return(out)
}

#' aggregate excursion tests spatially
#'
#' @param map.grid.cell from distanceGrid()
#' @param events an events object
#' @param agg.method What method do you want to use to aggregate the p-values? Choose from robustNull (default), fisher, sidak or lancaster
#' @param min.pval min p value
#'
#' @return gridded significance results
#' @export
#'
  spatialSigTest <- function(map.grid.cell,
                             events,
                             agg.method="robustNull",
                             min.pval = 0.001,
                             distance.cutoff = 2000,
                             use.weights = TRUE){

    TSids <- map.grid.cell$dist %>%
      dplyr::filter(distance < distance.cutoff)  %>%
 #     mutate(weight = 1 / (distance * 10/distance.cutoff)) %>%
      dplyr::distinct()

    if(nrow(TSids) == 0){
      map.grid.cell$pvalAbove <- NA
      map.grid.cell$pvalBelow <- NA
      map.grid.cell$pvalBoth <- NA
      map.grid.cell$pvalEither <- NA
      return(map.grid.cell)
    }

    #get the data within range for this grid cell
    resultsNow <- filter(events,paleoData_TSid %in% TSids$TSid) %>%
      rename(TSid = paleoData_TSid) %>%
      left_join(TSids,by = "TSid") %>%
      mutate(weight = gaspari_cohn(distance,distance.cutoff))

    if(agg.method != "robustNull"){
      uniquesWeights <- unique(resultsNow$weight)
      for (ii in uniquesWeights){
        sameWeightsLoc <- which(resultsNow$weight == ii)
        sameWeights <- resultsNow[sameWeightsLoc,]
        if (nrow(sameWeights)>1){
          resultsNow[sameWeightsLoc[1],]$pvalue_either <- aggregation::sidak(sameWeights$pvalue_either)
          resultsNow[sameWeightsLoc[-1],]$pvalue_either <- NA
          resultsNow[sameWeightsLoc[1],]$pvalue_both <- aggregation::sidak(sameWeights$pvalue_both)
          resultsNow[sameWeightsLoc[-1],]$pvalue_both <- NA
          resultsNow[sameWeightsLoc[1],]$pvalue_above <- aggregation::sidak(sameWeights$pvalue_above)
          resultsNow[sameWeightsLoc[-1],]$pvalue_above <- NA
          resultsNow[sameWeightsLoc[1],]$pvalue_below <- aggregation::sidak(sameWeights$pvalue_below)
          resultsNow[sameWeightsLoc[-1],]$pvalue_below <- NA
        }
      }
      resultsNow <- filter(resultsNow,!is.na(pvalue_either))
    }else{
      #only include the data with the lowest p-values for each site
      resultsNow$smallest <- FALSE
      uniquesWeights <- unique(resultsNow$weight)
      for (ii in uniquesWeights){
        sameWeightsLoc <- which(resultsNow$weight == ii)
        sameWeights <- resultsNow[sameWeightsLoc,]
        if (nrow(sameWeights)>1){
          smallest <- which.min(sameWeights$pvalue_either)
          resultsNow$smallest[sameWeightsLoc[smallest]] <- TRUE
        }else{
          resultsNow$smallest[sameWeightsLoc] <- TRUE
        }
      }
      resultsNow <- filter(resultsNow,smallest == TRUE)
    }

    if(min.pval > 0){
      resultsNow <- resultsNow %>%
        rowwise() %>%
        mutate(across(starts_with("pvalue_"),\(x) max(x,min.pval)))
    }

    if(nrow(resultsNow)==0){
      map.grid.cell$pval <- NA
    }else{
      if (agg.method=="fisher"){
        map.grid.cell$pvalAbove <- aggregation::fisher(resultsNow$pvalue_above)
        map.grid.cell$pvalBelow <- aggregation::fisher(resultsNow$pvalue_below)
        map.grid.cell$pvalBoth <- aggregation::fisher(resultsNow$pvalue_both)
        map.grid.cell$pvalEither <- aggregation::fisher(resultsNow$pvalue_either)

      }else if (agg.method=="sidak"){
        map.grid.cell$pvalAbove <- aggregation::sidak(resultsNow$pvalue_above)
        map.grid.cell$pvalBelow <- aggregation::sidak(resultsNow$pvalue_below)
        map.grid.cell$pvalBoth <- aggregation::sidak(resultsNow$pvalue_both)
        map.grid.cell$pvalEither <- aggregation::sidak(resultsNow$pvalue_either)
        }else if (agg.method=="lancaster"){
          map.grid.cell$pvalAbove <- aggregation::lancaster(resultsNow$pvalue_above,weights = resultsNow$weight)
          map.grid.cell$pvalBelow <- aggregation::lancaster(resultsNow$pvalue_below,weights = resultsNow$weight)
          map.grid.cell$pvalBoth <- aggregation::lancaster(resultsNow$pvalue_both,weights = resultsNow$weight)
          map.grid.cell$pvalEither <- aggregation::lancaster(resultsNow$pvalue_either,weights = resultsNow$weight)
      }else if (agg.method=="robustNull"){
        if(use.weights){
          resultsNow$weight <- map_dbl(resultsNow$weight,\(x) max(0.0001,min(x,1)))
          map.grid.cell <- append(map.grid.cell,calculateMultiTestSignificance(resultsNow,weights = resultsNow$weight))
        }else{
          map.grid.cell <- append(map.grid.cell,calculateMultiTestSignificance(resultsNow))
        }
      }else{
        stop("agg.method must be one of: fisher, sidak, lancaster or robustNull")
      }

    }
    map.grid.cell$resultsUsed <- resultsNow
    if(nrow(resultsNow) > 0){
      map.grid.cell$sumWeight <- sum(resultsNow$weight)
      map.grid.cell$maxWeight <- max(resultsNow$weight)
      top2 <- c(sort(resultsNow$weight,decreasing = TRUE),0)[1:2]

      map.grid.cell$sumTop2Weight <- sum(top2)

    }else{
      map.grid.cell$sumWeight <- NA
      map.grid.cell$maxWeight <- NA
      map.grid.cell$sumTop2Weight <- NA

    }

    return(map.grid.cell)
  }



#' plot results of spatialSigTest()
#'
#' @param pval.grid from spatialSigTest()
#' @param sigTestResults from excursionTestHighRes()
#' @param color.breaks at what significance levels should we put the color breaks?
#' @import ggplot2
#'
#' @return plot and plot data
#' @export
#'
plotSignificance <- function(pval.grid=NULL,
                             sigTestResults,
                             color.breaks = c(.001,.01,.05,.1,.2),
                             which.test = "pval",
                             restrict.sites = TRUE,
                             alpha.by.weight = TRUE,
                             cutoff.distance = 1500,
                             x.lim = c(-180,180),
                             y.lim = c(-90,90),
                             projection = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"){

  getPvals <- function(grid,which.test){
    if(is.null(grid[[which.test]])){
      return(NA)
    }else{
      return(grid[[which.test]])
    }
  }

  pvals <- map_dbl(pval.grid,getPvals,which.test)
  lat <- map_dbl(pval.grid,pluck,"lat")
  lon <- map_dbl(pval.grid,pluck,"lon")
  good <- which(!is.na(pvals))

  weight <- map_dbl(pval.grid[good],pluck,"sumTop2Weight")

  weight[weight < 0] <- 0



  plotData <- data.frame(pval = pvals[good],lat = lat[good],lon = lon[good],weight = weight)

  pvalOptions <- c("pvalPos","pvalNeg","pvalEither","pvalNet")


  handleNet <- function(positive,negative){
    if(positive > 0.5 & negative > 0.5){
      return(0.5)
    }

    out <- ifelse(positive > negative, 1 - negative,  positive)

    return(out)
  }

  sigTestResults$pvalPlot <- dplyr::case_when(
    which.test == "pvalPos" ~  sigTestResults$pvalue_positive,
    which.test == "pvalNeg" ~  sigTestResults$pvalue_negative,
    which.test == "pvalEither" ~  sigTestResults$pvalue_either,
    which.test == "pvalBoth" ~  sigTestResults$pvalue_both,
    which.test == "pvalNet" ~  purrr::map2_dbl(sigTestResults$pvalue_positive, sigTestResults$pvalue_negative, handleNet)
  )

  #optionally only plot the most significant
  if(restrict.sites){
    sigTestResults$plot <- FALSE

    udsn <- unique(sigTestResults$datasetId)
    for(u in udsn){
      w <- which(sigTestResults$datasetId == u)
      if(length(w) == 1){
        sigTestResults$plot[w] <- TRUE
      }else{
        if(which.test != "pvalNet"){
          tp <- w[which.min(sigTestResults$pvalPlot[w])]
          sigTestResults$plot[tp] <- TRUE
        }else{
          tp <- w[which.max(abs(sigTestResults$pvalPlot[w]-0.5))]
          sigTestResults$plot[tp] <- TRUE
        }
      }

    }
    sigTestResults <- filter(sigTestResults,plot == TRUE)

  }



  #sort by significance
  if(which.test == "pvalNet"){
  sigTestResults <- sigTestResults %>%
    arrange(abs(pvalPlot - 0.5))
  }else{
    sigTestResults <- sigTestResults %>%
      arrange(desc(pvalPlot))
  }

  cut <- 1.5
  world <- map_data("world") %>%
    filter(abs(long) < 180-cut,
           abs(lat) < 90-cut)

  toPlot <- filter(plotData,
                   !is.na(pval),
                   abs(lat) < 90-cut,
                   abs(lon) < 180-cut)


  p1 <- ggplot(data = toPlot,mapping = aes(x = lon,y = lat, fill = pval))


  if(alpha.by.weight){
    p1 <- p1 + geom_tile(aes(alpha = weight))
  }else{
    p1 <- p1 + geom_tile()
  }

  p1 <- p1 +
    geom_polygon(data=world, aes(x = long, y = lat, group = group), color="gray50", fill=alpha("white",0.1), inherit.aes = FALSE) +
    geom_point(inherit.aes = FALSE, data=sigTestResults, aes(x=geo_longitude, y=geo_latitude, fill=pvalPlot), size = 2, shape=21, color="black")+
    scale_fill_viridis_b(breaks=color.breaks)+
    scale_x_continuous("",breaks = seq(-180,180,by = 60)) +
    scale_y_continuous("",breaks = c(-87.5,seq(-60,60,by = 30),87.5)) +
    theme(legend.title = element_blank())+
    cowplot::theme_minimal_grid()+
    coord_sf(xlim=x.lim,
             ylim = y.lim,
             crs = projection,
             expand = TRUE,
             datum = sf::st_crs(4326),
             default_crs = sf::st_crs(4326))

  #print(p1)

  returns <- list(plot = p1, plot.data=plotData)

  return(returns)
}









