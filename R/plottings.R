actR_ggtheme <- ggplot2::theme_bw



#' Plot an excursion with uncertainties
#' @importFrom glue glue
#' @importFrom dplyr filter
#' @importFrom geoChronR plotTimeseriesEnsRibbons
#' @import ggplot2
#'
#' @param x The output of detectExcursion()
#' @param alpha What significance level to use?
#' @param print.significance Show significance on the plot?
#' @param x.axis.label Label the x-axis (default = NA, which will automatically generate from input)
#' @param y.axis.label Label the y-axis (default = NA, which will automatically generate from input)
#' @param lab.mult When labeling significance, how much higher than the highest point to put the label
#'
#' @return a ggplot2 object
plot.excursion <- function(x,
                           alpha = 0.05,
                           print.significance = FALSE,
                           x.axis.label = NA,
                           y.axis.label = NA,
                           lab.mult = 0.02){

  if(is.na(x.axis.label)){
    x.axis.label <- glue::glue("{x$timeVariableName} ({x$timeUnits})")
  }
  if(is.na(y.axis.label)){
    y.axis.label <- glue::glue("Detrended {x$paleoData_variableName} ({x$paleoData_units})")
  }

  #determine what is ensemble
  ensOut <- x$eventDetection[[1]]
  hasTimeEnsemble <- !identicalVectorsList(ensOut$time)
  hasPaleoEnsemble <- !identicalVectorsList(ensOut$vals)

  if(hasTimeEnsemble){
    timeMat <- list2matrix(ensOut$time)
  }else{
    timeMat <- ensOut$time[[1]]
  }

  if(hasPaleoEnsemble){
    paleoMat <- list2matrix(ensOut$vals)
  }else{
    paleoMat <- ensOut$vals[[1]]
  }

  #if time or values are ensemble, plot as ribbons, with single line over top.
  if(hasPaleoEnsemble | hasTimeEnsemble){
    ribbons <- geoChronR::plotTimeseriesEnsRibbons(X = timeMat,Y = paleoMat,probs = c(.025,.25,.75,.975))

    #check for significance
    if(x$empirical_pvalue < alpha){
      alpha.msg <- glue::glue("Empirical p-value {x$empirical_pvalue} is < {alpha}")
      ts <- dplyr::filter(ensOut,eventDetected == TRUE)
    }else{
      alpha.msg <- glue::glue("Empirical p-value {x$empirical_pvalue} exceeds {alpha}")
      ts <- dplyr::filter(ensOut,eventDetected == FALSE)
    }

    #pick a representative ensemble member
    rm <- sample(which(ts$nExcusionVals == median(ts$nExcusionVals,na.rm = TRUE)),size = 1)
    em <- ts[rm,]

    plotout <- plot.excursionCore(em,add.to.plot = ribbons)

    if(print.significance){
      l.y <- max(em$vals[[1]],na.rm = TRUE)+lab.mult*abs(diff(range(em$vals[[1]]),na.rm = TRUE))
      l.x <- mean(em$time[[1]],na.rm = TRUE)
      plotout <- plotout +
        ggplot2::annotate("text",x = l.x, y = l.y,label = alpha.msg)
    }
  }else{
    plotout <- plotExcursionCore(ensOut[1,])
  }


  #title
  if(!is.na(x$dataSetName) & !is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$dataSetName} - {x$paleoData_variableName}: Excursion")
  }else if(is.na(x$dataSetName) & !is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$paleoData_variableName}: Excursion")
  }else if(!is.na(x$dataSetName) & is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$dataSetName}: Excursion")
  }else{
    title <- glue::glue("Excursion")
  }

  #add labels, directionality
  plotout <- plotout +
    xlab(x.axis.label)+
    ylab(y.axis.label)+
    ggtitle(title)


  if(grepl(x = x$timeUnits,, pattern = "ky",ignore.case = TRUE) | grepl(x = x$timeUnits, pattern = "bp",ignore.case = TRUE)){
    plotout <- plotout +
      scale_x_reverse()
  }

  return(plotout)
}



#' Plot a basic excursion without propagated uncertainty
#'
#' @import ggplot2
#' @param x the output of detectExcursionCore()
#' @param add.to.plot a ggplot object upon which to add this plot
#'
#' @return a ggplot object
plot.excursionCore <- function(x,
                               add.to.plot = ggplot2::ggplot()){

  #spread out the list columns
  tp <- tidyr::unchop(x, c("time","vals","isExcursion"))

  params <- createTibbleFromParameterString(x$parameters[1])

  #Prepare some columns for plotting
  pre_i <- which(tp$time < tp$time_start)
  post_i <- which(tp$time > tp$time_end)

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
    geom_point(data = tp,aes(x = time, y = vals,color = isExcursion)) +
    geom_line(data = tp,aes(x = time, y = vals)) +
    geom_line(data = tp,aes(x = time, y = refMeans), color = 'red',linetype = 1)+
    geom_line(data = tp,aes(x = time, y = refSdHi),  color = 'red', linetype = 4)+
    geom_line(data = tp,aes(x = time, y = refSdLo),  color = 'red', linetype = 4)+
    ylab("Detrended values")+
    scale_color_manual("Excursion Detected",values = c("black","red"))+
    actR_ggtheme()

  return(plotOut)
}

#' Plot a shift in mean or variance
#' @import tibble ggplot2 dplyr
#'
#' @param x Output of detectShiftCore()
#' @param time Input time vector
#' @param vals input value vector
#' @param add.to.plot ggplot to add this to
#' @param mean.color color of the mean lines
#'
#' @return A ggplot object
plotSectionMeans <- function(x,time,vals,add.to.plot = ggplot2::ggplot(),mean.color = "red"){

  #find section means
  cpt <- sort(c(min(time,na.rm = TRUE),x$time_start,max(time,na.rm = TRUE)))

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

#' Plot a shift in mean or variance
#' @import tibble ggplot2 dplyr
#' @param line.color color of the line
#' @param mean.color color of the means
#' @param x tibble with time and vals
#'
#' @return A ggplot object
plot.shiftCore <- function(x,line.color = "black", mean.color = "red"){
  if(nrow(x)==0){#how to handle this?

  }
  time <- x$time[[1]]
  vals <- x$vals[[1]]

  #find section means
  cpt <- sort(c(min(time,na.rm = TRUE),x$time_start,max(time,na.rm = TRUE)))

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

#' Plot mean or variance shifts, with uncertainties and null hypothesis testing
#'
#' @importFrom geoChronR plotTimeseriesEnsRibbons
#' @import ggplot2 tidyr RColorBrewer purrr dplyr egg
#'
#' @param x Output from actR::detectShift
#' @param cl.color Color palette, single color, or vector of colors to use for confidence intervals (default = "Reds"s)
#' @param plot.sig.vlines Plot vertical lines at significant change points? (default = TRUE)
#' @param label.sig Label significant change points (default = TRUE)
#' @param alpha significance level (default = 0.05)
#' @param x.axis.label Label the x-axis (default = NA, which will automatically generate from input)
#' @param y.axis.label Label the y-axis (default = NA, which will automatically generate from input)
#' @param combine.plots Combine the probability and timeseries plots into a single plot (TRUE)? Or return a list with each plot as a separate object (FALSE)?
#' @param x.lims 2-element vector to usee as x-axis limits (default = NA)
#'
#' @return a ggplot object
plot.shift <- function(x,
                       cl.color = "Reds",
                       plot.sig.vlines = TRUE,
                       label.sig = TRUE,
                       alpha = 0.05,
                       x.axis.label = NA,
                       x.lims = NA,
                       y.axis.label = NA,
                       combine.plots = TRUE){

  if(is.na(x.axis.label)){
    x.axis.label <- glue::glue("{x$input$timeVariableName} ({x$input$timeUnits})")
  }
  if(is.na(y.axis.label)){
    y.axis.label <- glue::glue("{x$input$paleoData_variableName} ({x$input$paleoData_units})")
  }


  paramTib <- createTibbleFromParameterString(x$parameters[1])

  #plot ensemble ribbons
  ribbons <- geoChronR::plotTimeseriesEnsRibbons(X = x$timeEns,Y = x$valEns) + actR_ggtheme()

  #get shift type
  if(grepl(pattern = "cpt.mean",paramTib$cpt.fun,ignore.case = T) & !grepl(pattern = "cpt.meanVar",paramTib$cpt.fun,ignore.case = T)){
    shift.type <- "Shift in Mean"
  }else if(grepl(pattern = "cpt.var",paramTib$cpt.fun,ignore.case = T)){
    shift.type <- "Shift in Variance"
  }else if(grepl(pattern = "cpt.meanVar",paramTib$cpt.fun,ignore.case = T)){
    shift.type <- "Shift in Mean and Variance"
  }else{
    shift.type <- paramTib$cpt.fun
  }
  #title
  if(!is.na(x$input$dataSetName) & !is.na(x$input$paleoData_variableName)){
    title <- glue::glue("{x$input$dataSetName} - {x$input$paleoData_variableName}: {shift.type}")
  }else if(is.na(x$input$dataSetName) & !is.na(x$input$paleoData_variableName)){
    title <- glue::glue("{x$input$paleoData_variableName}: {shift.type}")
  }else if(!is.na(x$input$dataSetName) & is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$input$dataSetName}: {shift.type}")
  }else{
    title <- glue::glue("{shift.type}")
  }

  #add title and labels



  timeMed <- apply(x$timeEns,1,median,na.rm = TRUE)
  valMed <- apply(x$valEns,1,median,na.rm = TRUE)

  #plot changepoint probabilities
  #make step plot

  cpp <- x$shiftDetection %>%
    dplyr::mutate(time_end = time_end - .001) %>%
    tidyr::pivot_longer(c("time_start","time_end"),values_to = "time_edges") %>%
    dplyr::arrange(time_edges)

  #get x.range
  if(is.na(x.lims)){
    x.lims <- range(cpp$time_edges)
  }

  npp <- cpp %>%
    tidyr::pivot_longer(starts_with("q"),names_to = "cl",values_to = "nullProbs")

  np <- length(unique(npp$cl))
  #deal with line colors
  if(cl.color %in% rownames(RColorBrewer::brewer.pal.info)){#
    #then it's an RColorBrewer pallette
    colorScale <- rep_len(suppressWarnings(RColorBrewer::brewer.pal(n = np,name = cl.color)),length.out = np)
  }else{#it's not
    if(length(cl.color) == 1){#apply one color to all
      colorScale <- rep(cl.color,times = np)
    }else{
      if(length(cl.color) == np){
        colorScale <- cl.color
      }else{
        stop("cl.color must be either 1) a single color to repeated, 2) an RColorBrewer palette or 3) a string the same length of the number of lines to be plotted")
      }
    }
  }

  #get significant events
  sig.event <- summarizeShiftSignificance(x$shiftDetection, alpha = alpha, paramTib = paramTib )

  if(nrow(sig.event) == 0){
    any.sig = FALSE
  }else{
    any.sig = TRUE
  }

  minp <- x$shiftDetection %>%
    dplyr::filter(empirical_pvalue > 0) %>%
    dplyr::select("empirical_pvalue") %>%
    min(na.rm = TRUE)

  max.y <- x$shiftDetection %>%
    dplyr::select("event_probability" | starts_with("q")) %>%
    max(na.rm = TRUE)

  sig.event$pvallab <- paste("p =",sig.event$empirical_pvalue)
  sig.event$pvallab[sig.event$empirical_pvalue == 0] <- glue::glue("p < {minp}")

  #plot shift frequency and cls
  probPlot <- ggplot()+
    actR_ggtheme()+
    geom_ribbon(data = cpp,aes(x = time_edges,ymin = 0,ymax = event_probability)) +
    geom_line(data = npp,aes(x = time_edges,y = nullProbs,color = cl))+
    scale_color_manual(values = colorScale)+
    xlab(x.axis.label) +
    ylab("Shift Frequency")+
    theme(legend.position = c(0.9,0.8),
          legend.background = element_blank(),
          legend.title = element_blank())


  if(grepl(x = x$input$timeUnits,, pattern = "ky",ignore.case = TRUE) | grepl(x = x$input$timeUnits, pattern = "bp",ignore.case = TRUE)){
    this_x_scale <- scale_x_reverse
    x.lims <- rev(x.lims)
  }else{
    this_x_scale <- scale_x_continuous
  }

  probPlot <- probPlot +
    this_x_scale(limits = x.lims)


  if(plot.sig.vlines & any.sig){
    probPlot <- probPlot +
      geom_vline(data = sig.event,aes(xintercept = time_mid),linetype = "dashed",color = "gray50")
  }
  if(label.sig & any.sig){
    probPlot <- probPlot +
      geom_label(data = sig.event,aes(x = time_mid,y = max.y,label = pvallab))
  }


  timeSeries <- plotSectionMeans(add.to.plot = ribbons,sig.event,time = timeMed,vals = valMed)+
    this_x_scale(name = x.axis.label, position = "top")+
    scale_y_continuous(name = y.axis.label, position = "right")+
    theme(axis.ticks.x.bottom = element_blank(),
          plot.margin=unit(c(1,1,-.2,1), "cm"))+
    ggtitle(title)+
    coord_cartesian(xlim =  x.lims)



  if(combine.plots){
    # combine the two plots
    out <- egg::ggarrange(plots = list(timeSeries,probPlot),ncol = 1,padding = 0)
  }else{
    out <- list(timeSeriesPlot = timeSeries, probabilityPlot = probPlot)
  }
  return(out)

}



