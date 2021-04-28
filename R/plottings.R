actR_ggtheme <- ggplot2::theme_bw



#' Plot an excursion with uncertainties
#'
#' @param exc.out The output of detectExcursion()
#'
#' @return a ggplot2 object
#' @export
plotExcursion <- function(exc.out){

  #determine what is ensemble
  ensOut <- exc.out$eventDetection[[1]]
  hasAgeEnsemble <- !identicalVectorsList(ensOut$age)
  hasPaleoEnsemble <- !identicalVectorsList(ensOut$vals)

  if(hasAgeEnsemble){
    ageMat <- list2matrix(ensOut$age)
  }else{
    ageMat <- ensOut$age[[1]]
  }

  if(hasPaleoEnsemble){
    paleoMat <- list2matrix(ensOut$vals)
  }else{
    paleoMat <- ensOut$vals[[1]]
  }

  #if age or values are ensemble, plot as ribbons, with single line over top.
  if(hasPaleoEnsemble | hasAgeEnsemble){
    ribbons <- geoChronR::plotTimeseriesEnsRibbons(X = ageMat,Y = paleoMat,probs = c(.025,.25,.75,.975))

    #pick a representative ensemble member
    rm <- sample(which(ensOut$nExcusionVals == median(ensOut$nExcusionVals,na.rm = TRUE)),size = 1)
    em <- ensOut[rm,]

    plotout <- plotExcursionCore(em,add.to.plot = ribbons)
  }else{
    plotout <- plotExcursionCore(ensOut[1,])
  }

  return(plotout)
}



#' Plot a basic excursion without propagated uncertainty
#'
#' @param exc.out.core the output of detectExcursionCore()
#' @param add.to.plot a ggplot object upon which to add this plot
#'
#' @return a ggplot object
#' @export
plotExcursionCore <- function(exc.out.core,
                              add.to.plot = ggplot2::ggplot()){

  #spread out the list columns
  tp <- tidyr::unchop(exc.out.core, c("age","vals","isExcursion"))

  params <- createTibbleFromParameterString(exc.out.core$parameters[1])

  #Prepare some columns for plotting
  pre_i <- which(tp$age < tp$time_start)
  post_i <- which(tp$age > tp$time_end)

  tp$refMeans <- NA
  tp$refMeans[pre_i] <- tp$preMean[pre_i]
  tp$refMeans[post_i] <- tp$postMean[post_i]

  tp$refSdHi <- NA
  tp$refSdHi[pre_i] <- tp$refMeans[pre_i] + tp$preSd[pre_i] * params$sig.num
  tp$refSdHi[post_i] <- tp$refMeans[post_i] + tp$postSd[post_i] * params$sig.num

  tp$refSdLo <- NA
  tp$refSdLo[pre_i] <- tp$refMeans[pre_i] - tp$preSd[pre_i] * params$sig.num
  tp$refSdLo[post_i] <- tp$refMeans[post_i] - tp$postSd[post_i] * params$sig.num

  #Find the excursions
  plotOut <- add.to.plot +
    geom_point(data = tp,aes(x = age, y = vals,color = isExcursion)) +
    geom_line(data = tp,aes(x = age, y = vals)) +
    geom_line(data = tp,aes(x = age, y = refMeans), color = 'red',linetype = 1)+
    geom_line(data = tp,aes(x = age, y = refSdHi),  color = 'red', linetype = 4)+
    geom_line(data = tp,aes(x = age, y = refSdLo),  color = 'red', linetype = 4)+
    ylab("Detrended values")+
    scale_color_manual("Excursion Detected",values = c("black","red"))+
    actR_ggtheme()

  return(plotOut)
}

plotSectionMeans <- function(ms.core.out,time,vals,add.to.plot = ggplot2::ggplot()){

  #find section means
  cpt <- sort(c(min(time,na.rm = TRUE),ms.core.out$time_start,max(time,na.rm = TRUE)))

  ind <- means <- matrix(NA,nrow = length(time))

  for(i in 1:(length(cpt)-1)){
    w <- which(dplyr::between(time,cpt[i],cpt[i+1]))
    ind[w] <- i
    means[w] <- mean(vals[w],na.rm = TRUE)
  }

  tp <- tibble::tibble(time,
                       vals,
                       cpt_section = factor(ind),
                       cpt_mean = means)


  #plot it. This is kinda clumsy for the moment
  ms <- add.to.plot+
    actR_ggtheme()


  for(i in 1:(length(cpt)-1)){
    td <- dplyr::filter(tp,cpt_section == i)
    ms <- ms + geom_line(data = td,aes(time,cpt_mean),color = mean.color)
  }

  return(ms)
}


plotMeanShiftCore <- function(ms.core.out,line.color = "black", mean.color = "red"){
  if(nrow(ms.core.out)==0){#how to handle this?

  }
  time <- ms.core.out$age[[1]]
  vals <- ms.core.out$vals[[1]]

#find section means
  cpt <- sort(c(min(time,na.rm = TRUE),ms.core.out$time_start,max(time,na.rm = TRUE)))

  ind <- means <- matrix(NA,nrow = length(time))

  for(i in 1:(length(cpt)-1)){
    w <- which(dplyr::between(time,cpt[i],cpt[i+1]))
    ind[w] <- i
    means[w] <- mean(vals[w],na.rm = TRUE)
  }

tp <- tibble::tibble(time,
                     vals,
                     cpt_section = factor(ind),
                     cpt_mean = means)


#plot it. This is kinda clumsy for the moment
ms <- ggplot()+
  geom_line(data = tp,aes(time,vals),color = line.color)+
  actR_ggtheme()


for(i in 1:(length(cpt)-1)){
  td <- dplyr::filter(tp,cpt_section == i)
  ms <- ms + geom_line(data = td,aes(time,cpt_mean),color = mean.color)
}

return(ms)

}

plotMeanShift <- function(ms.out){

  paramTib <- createTibbleFromParameterString(ms.out$parameters)

  #plot ensemble ribbons
  ribbons <- geoChronR::plotTimeseriesEnsRibbons(X = ms.out$ageEns,Y = ms.out$valEns) + actR_ggtheme()

  ageMed <- apply(ms.out$ageEns,1,median,na.rm = TRUE)
  valMed <- apply(ms.out$valEns,1,median,na.rm = TRUE)

  #plot changepoint probabilities
  #make step plot
  cpp <- ms.out$shiftDetection %>%
    tidyr::pivot_longer(c("time_start","time_end"),values_to = "time_edges")

  npp <- cpp %>%
    tidyr::pivot_longer(starts_with("q"),names_to = "cl",values_to = "nullProbs")

  probPlot <- ggplot()+
    actR_ggtheme()+
    geom_ribbon(data = cpp,aes(x = time_edges,ymin = 0,ymax = event_probability)) +
    geom_line(data = npp,aes(x = time_edges,y = nullProbs,color = cl))+
    scale_color_brewer(palette = "Set1")+
    ylab("Shift Frequency")+
    theme(legend.position = c(0.9,0.8),
          legend.background = element_blank(),
          legend.title = element_blank())



  ms.out$shiftDetection


  #find significant events
    sig.event <- ms.out$shiftDetection %>%
      dplyr::filter(event_probability > `q0.95`) %>%
      dplyr::mutate(exceedance = event_probability - `q0.95`) %>%
      dplyr::arrange(dplyr::desc(exceedance))

    #remove those within minimum distance
    bm <- sig.event$time_mid
    tr <- c()
    for(i in 1:(length(bm) - 1)){
      diffs <- abs(bm-bm[i])
      close <- which(diffs < paramTib$minimum.segment.length)
      ttr <- intersect(close,seq(i+1,length(bm)))
      tr <- c(tr,ttr)
    }

    if(length(tr) > 0){
      sig.event <- sig.event[-tr,]
    }

    timeSeries <- plotSectionMeans(add.to.plot = ribbons,sig.event,time = ageMed,vals = valMed)+
      scale_x_continuous(position = "top")+scale_y_continuous(position = "right")+
     theme(axis.ticks.x.bottom = element_blank(),
            plot.margin=unit(c(1,1,-.2,1), "cm"))

  # combine the two plots
    out <- egg::ggarrange(plots = list(timeSeries,probPlot),ncol = 1,padding = 0)

    return(out)

}



