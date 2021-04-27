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

plotMeanShift <- function(ms.out,prob.thresh = .5){
#plot ensemble ribbons
  #determine what is ensemble
  ensOut <- ms.out$eventDetection[[1]]
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
  #plot change point probability

#two options
  #1. find an ensemble member with change points most similar to the median probabilities
  #2. calculate means from ensemble, and show along ensemblee median.


}
