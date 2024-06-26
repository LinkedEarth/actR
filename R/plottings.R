#' @export
#' @family plot help
#' @title Define a plot theme for actR
#' @description Use this to define a theme across actR
#' @import ggplot2
#' @param font.family Specify a font family to use for the theme (default = "Helvetica")
#' @param ... parameters to pass to theme function
actR_ggtheme <- function(font.family = "Helvetica",...){
  ggplot2::theme_bw(base_family = font.family,...)
}


#' Plot an excursion with uncertainties
#'
#' @param x The output of detectExcursion()
#' @param ... additional parameters (see plotExcursion)
#'
#' @return a ggplot2 object
#' @export
plot.excursion <- function(x,...){
  return(plotExcursion(x,...))
}

#' Plot an excursion with uncertainties
#' @importFrom glue glue
#' @importFrom dplyr filter
#' @importFrom geoChronR plotTimeseriesEnsRibbons
#' @import ggplot2
#' @param x The output of detectExcursion()
#' @param alpha What significance level to use?
#' @param print.significance Show significance on the plot?
#' @param x.axis.label Label the x-axis (default = NA, which will automatically generate from input)
#' @param y.axis.label Label the y-axis (default = NA, which will automatically generate from input)
#' @param lab.mult When labeling significance, how much higher than the highest point to put the label
#'
#' @return a ggplot2 object
#' @export
plotExcursion <- function(x,
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
  ensOut <- x$event_detection[[1]]
  hasTimeEnsemble <- suppressWarnings(!identicalVectorsList(ensOut$time))
  hasPaleoEnsemble <- suppressWarnings(!identicalVectorsList(ensOut$vals))
  hasParamEnsemble <- hasParameterEnsemble(x)

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

  if(hasParamEnsemble){
    if(hasPaleoEnsemble){
      #limit outliers
      m95 <- quantile(paleoMat,na.rm = TRUE,c(.025,.975))
      bad <- which(paleoMat < m95[1] | paleoMat > m95[2])
      pmb <- paleoMat
      pmb[bad] <- NA
      nGoodInCol <- colSums(!is.na(pmb))
      badCol <- which(nGoodInCol < quantile(nGoodInCol,probs = .25))
      paleoMat <- paleoMat[,-badCol]
    }
    if(hasTimeEnsemble){
      timeMat <- timeMat[,-badCol]
    }
  }


  #if time or values are ensemble, plot as ribbons, with single line over top.
  if(hasPaleoEnsemble | hasTimeEnsemble){
    ribbons <- geoChronR::plotTimeseriesEnsRibbons(X = timeMat,
                                                   Y = paleoMat,
                                                   probs = c(.025,.25,.75,.975))

    #check for significance
    if(x$pvalue < alpha){
      alpha.msg <- glue::glue("Empirical p-value {x$pvalue} is < {alpha}")
      ts <- dplyr::filter(ensOut,eventDetected == TRUE)
    }else{
      alpha.msg <- glue::glue("Empirical p-value {x$pvalue} exceeds {alpha}")
      ts <- dplyr::filter(ensOut,eventDetected == FALSE)
    }

    if(nrow(ts) == 0){
      stop("No representative ensemble member found, which is a problem. It probably means you should increase n.ens and maybe null.hypothesis.n in detectExcursion.")
    }

    #pick a representative ensemble member
    sami <- which(ts$nExcursionVals == median(ts$nExcursionVals,na.rm = TRUE))
    if(length(sami) == 0){
      sami <- seq_along(ts$nExcursionVals)
    }
    rm <- sample(sami,size = 1)

    em <- ts[rm,]

    plotout <- plot.excursionCore(em,add.to.plot = ribbons)

    if(print.significance){
      l.y <- max(em$vals[[1]],na.rm = TRUE)+lab.mult*abs(diff(range(em$vals[[1]]),na.rm = TRUE))
      l.x <- mean(em$time[[1]],na.rm = TRUE)
      plotout <- plotout +
        ggplot2::annotate("text",x = l.x, y = l.y,label = alpha.msg)
    }
  }else{
    plotout <- plot.excursionCore(ensOut[1,])
  }


  #title
  if(!is.na(x$dataSetName) & !is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$dataSetName} - {x$paleoData_variableName}: Excursion ({x$exc.type[1]})")
  }else if(is.na(x$dataSetName) & !is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$paleoData_variableName}: Excursion ({x$exc.type[1]})")
  }else if(!is.na(x$dataSetName) & is.na(x$paleoData_variableName)){
    title <- glue::glue("{x$dataSetName}: Excursion ({x$exc.type[1]})")
  }else{
    title <- glue::glue("Excursion ({x$exc.type[1]})")
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
#' @param ... a ggplot object upon which to add this plot
#'
#' @return a ggplot object
#' @export
plot.excursionCore <- function(x,...){
  return(plotExcursionCore(x,...))
}



#' Plot a basic excursion without propagated uncertainty
#'
#' @import ggplot2
#' @param x the output of detectExcursionCore()
#' @param add.to.plot a ggplot object upon which to add this plot
#' @importFrom tidyr unchop
#' @return a ggplot object
#' @export
plotExcursionCore <- function(x,
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
    geom_point(data = tp,aes(x = time, y = vals,color = isExcursion,shape = isExcursion)) +
    geom_line(data = tp,aes(x = time, y = vals)) +
    geom_line(data = tp,aes(x = time, y = refMeans), color = 'red',linetype = 1)+
    geom_line(data = tp,aes(x = time, y = refSdHi),  color = 'red', linetype = 4)+
    geom_line(data = tp,aes(x = time, y = refSdLo),  color = 'red', linetype = 4)+
    ylab("Detrended values")+
    scale_color_manual("Excursion Detected",values = c("black","red"))+
    scale_shape_manual("Excursion Detected",values = c(15,16))+
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
#' @export
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

#' Plot mean or variance shifts, with uncertainties and null hypothesis testing
#'
#' @param x Output from actR::detectShift
#' @param ... additional inputs (see plotShiftCore)
#' @export
#' @return A ggplot object
plot.shiftCore <- function(x,...){
  return(plotShiftCore(x,...))
}

#' Plot a shift in mean or variance
#' @import tibble ggplot2 dplyr
#' @param line.color color of the line
#' @param mean.color color of the means
#' @param x tibble with time and vals
#' @export
#' @return A ggplot object
plotShiftCore <- function(x,line.color = "black", mean.color = "red"){
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
#' @param x Output from actR::detectShift
#' @param ... more inputs, see (plotShift)
#'
#' @export
#' @return a ggplot object
plot.shift <- function(x,...){
  return(plotShift(x,...))
}


#' Plot mean or variance shifts, with uncertainties and null hypothesis testing
#'
#' @importFrom geoChronR plotTimeseriesEnsRibbons
#' @import ggplot2  RColorBrewer purrr dplyr egg
#' @importFrom tidyr pivot_longer
#'
#' @param x Output from actR::detectShift (tibble)
#' @param cl.color Color palette, single color, or vector of colors to use for confidence intervals (default = "Reds"s)
#' @param sf.fill Color palette, single color, or vector of colors to use for confidence intervals (default = c("grey"))
#' @param plot.sig.vlines Plot vertical lines at significant change points? (default = TRUE)
#' @param label.sig Label significant change points (default = TRUE)
#' @param qs significance levels to plot. The first value will be used to calculate the alpha for significant shifts (default = c(0.95,0.9))
#' @param x.axis.label Label the x-axis (default = NA, which will automatically generate from input)
#' @param y.axis.label Label the y-axis (default = NA, which will automatically generate from input)
#' @param x.lims 2-element vector to use as x-axis limits (default = NA)
#' @param y.lims.direction A value of -1 will flip the y-axis for both panels (default = 1) #TODO
#' @param shift.direction plot positive, negative or both ("positive/negative" - the default)
#' @param combine.plots Combine the probability and timeseries plots into a single plot (TRUE)? Or return a list with each plot as a separate object (FALSE)?
#'
#' @inheritDotParams geoChronR::plotTimeseriesEnsRibbons
#' @export
#' @return a ggplot object
plotShift <- function(x,
                       cl.color = "Reds",
                       sf.fill = c('navy','firebrick'), #todo
                       plot.sig.vlines = FALSE,
                       label.sig = FALSE,
                       qs = c(0.95,0.9),
                       x.axis.label = NA,
                       x.lims = NA,
                       y.axis.label = NA,
                       y.lims.direction = 1,
                       combine.plots = TRUE,
                       shift.direction = "positive/negative",
                       ...){



  if(any(is.na(x.axis.label))){
    x.axis.label <- glue::glue("{unique(x$time_variableName)} ({unique(x$timeUnits)})")
  }
  if(any(is.na(y.axis.label))){
    y.axis.label <- glue::glue("{unique(x$paleoData_variableName)} ({unique(x$paleoData_units)})")
  }


  #paramTib <-  (x$parameters[1])

  #plot ensemble ribbons
  ribbons <- geoChronR::plotTimeseriesEnsRibbons(X = x$time[[1]],Y = x$paleoData_values[[1]],...) + actR_ggtheme()

  #get shift type
  if(grepl(pattern = "cpt.mean",x$cpt.fun[[1]],ignore.case = T) & !grepl(pattern = "cpt.meanVar",x$cpt.fun[[1]],ignore.case = T)){
    shift.type <- "Shift in Mean"
  }else if(grepl(pattern = "cpt.var",x$cpt.fun[[1]],ignore.case = T)){
    shift.type <- "Shift in Variance"
  }else if(grepl(pattern = "cpt.meanVar",x$cpt.fun[[1]],ignore.case = T)){
    shift.type <- "Shift in Mean and Variance"
  }else{
    shift.type <- x$cpt.fun[[1]]
  }

  #title
  if(!any(is.na(x$dataSetName)) & !any(is.na(x$paleoData_variableName))){
    title <- glue::glue("{x$dataSetName} - {x$paleoData_variableName}: {shift.type}") |> unique()
  }else if(suppressWarnings(any(is.na(x$input$dataSetName)) & !any(is.na(x$paleoData_variableName)))){
    title <- glue::glue("{x$paleoData_variableName}: {shift.type}") |> unique()
  }else if(suppressWarnings(!any(is.na(x$input$dataSetName)) & any(is.na(x$paleoData_variableName)))){
    title <- glue::glue("{x$dataSetName}: {shift.type}") |> unique()
  }else{
    title <- glue::glue("{shift.type}")
  }

  #Median line
  timeMed <- apply(x$time[[1]],1,median,na.rm = TRUE)
  valMed <- apply(x$paleoData_values[[1]],1,median,na.rm = TRUE)

  #get x.range
  if(any(is.na(x.lims))){
    x.lims <- c(min(x$time_start),max(x$time_end))
  }

  #Create dataframe of probabilities to plot
  if (shift.direction == 'positive/negative'){
    cpp  <- x[,c("time_start","time_end","time_mid","event_probability_positive","null_probability_positive","pvalue_positive")]
    cpp2 <- x[,c("time_start","time_end","time_mid","event_probability_negative","null_probability_negative","pvalue_negative")]
  }else{
    cpp <- x[,c("time_start","time_end","time_mid",paste0("event_probability_",shift.direction),paste0("null_probability_",shift.direction),paste0("pvalue_",shift.direction))]
    cpp2 <- NULL
  }
  colnames(cpp) <- c("time_start","time_end","time_mid","event_probability","null_probability","pvalue")
  if(!is.null(cpp2)){colnames(cpp2) <- c("time_start","time_end","time_mid","event_probability","null_probability","pvalue")}

  #Add in null quantiles
  for (q in qs){
    cpp[[paste0('q',q)]] <- unlist(lapply(cpp$null_probability,quantile,q))
    if(!is.null(cpp2)){cpp2[[paste0('q',q)]] <- unlist(lapply(cpp2$null_probability,quantile,q))}
  }
  npp <- cpp %>%
    tidyr::pivot_longer(c("time_start","time_end"),values_to = "time_edges") %>%
    tidyr::pivot_longer(starts_with("q"),names_to = "cl",values_to = "nullProbs")
  if(!is.null(cpp2)){
    npp2 <- cpp2 %>%
      tidyr::pivot_longer(c("time_start","time_end"),values_to = "time_edges") %>%
      tidyr::pivot_longer(starts_with("q"),names_to = "cl",values_to = "nullProbs")
  }

  #deal with line colors
  np <- length(qs)
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

  #plot shift frequency and cls
  probPlot <- ggplot()+
    actR_ggtheme()+
    geom_hline(yintercept = 0)

  #Plot the data for the paleoData
  if (shift.direction == 'positive/negative'){
    probPlot <- probPlot +
      geom_col(data = cpp, aes(x=time_mid, y=event_probability,fill='positive'),
               width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=0.5)+
      #geom_col(data = cpp%>%filter(pvalue<1-qs[1]), aes(x=time_mid, y=event_probability,fill='positive'),
      #         width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=1)+
      geom_col(data = cpp[which(cpp$event_probability >= cpp[paste0('q',qs[1])]),], aes(x=time_mid, y=event_probability,fill='positive'),
               width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=1)+
      geom_col(data = cpp2, aes(x=time_mid, y=event_probability*-1,fill='negative'),
               width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=0.5)+
      geom_col(data = cpp2[which(cpp2$event_probability >= cpp2[paste0('q',qs[1])]),], aes(x=time_mid, y=-event_probability,fill='negative'),
               width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=1)+
      geom_path(data = npp,aes(x = time_edges, y = nullProbs, color = cl)) +
      geom_path(data = npp2,aes(x = time_edges, y = nullProbs*-1, color = cl))


  } else{
    probPlot <- probPlot +
      geom_col(data = cpp, aes(x=time_mid, y=event_probability,fill=shift.direction),
               width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=0.5)+
      geom_col(data = cpp[which(cpp$event_probability >= cpp[paste0('q',qs[1])]),], aes(x=time_mid, y=event_probability,fill=shift.direction),
              width=median(diff(x$time_mid)), color='black', linewidth=0.1, alpha=1)+
      geom_path(data = npp,aes(x = time_edges,y = nullProbs,color = cl))
  }


  probPlot <- probPlot +
    scale_color_manual(values = colorScale,name='Confidence\nLevel') +
    scale_fill_manual(values = sf.fill,name='Direction') #+
    #theme(#legend.position = c(0.9,0.8),
      #legend.background = element_rect(fill='white',color='black'),
      #legend.title = element_blank())

  #get significant events
  cpp_sig <- rbind(cpp,cpp2)
  cpp_sig$minimum.segment.length <- x$minimum.segment.length[1]
  sig.event <- summarizeShiftSignificance(cpp_sig, alpha = 1-qs[1],x$minimum.segment.length[1])

  if(nrow(sig.event) == 0){
    any.sig = FALSE
  }else{
    any.sig = TRUE
  }

  if(max(cpp_sig$pvalue,na.rm = TRUE) > 0){
    minp <- cpp_sig %>%
      dplyr::select(pvalue) %>%
      dplyr::filter(pvalue > 0) %>%
      min(na.rm = TRUE)
  }else{
    minp <- 0
  }

  max.y <- cpp %>%
    dplyr::select("event_probability" | starts_with("q")) %>%
    max(na.rm = TRUE)

  sig.event$pvallab <- paste("p =",sig.event$pvalue)
  sig.event$pvallab[sig.event$pvalue == 0] <- glue::glue("p < {minp}")

  if(grepl(x = x$timeUnits[[1]], pattern = "ky",ignore.case = TRUE) | grepl(x = x$timeUnits[[1]], pattern = "bp",ignore.case = TRUE)){
    this_x_scale <- scale_x_reverse
    x.lims <- rev(sort(x.lims))
  }else{
    this_x_scale <- scale_x_continuous
  }

  if(y.lims.direction<0){
    this_y_scale <- scale_y_reverse
  }else{
    this_y_scale <- scale_y_continuous
  }

  probPlot <- probPlot +
    this_y_scale(name="Shift Frequency")+
    this_x_scale(name=x.axis.label,limits = x.lims)


  if(plot.sig.vlines & any.sig){
    probPlot <- probPlot +
      geom_vline(data = sig.event,aes(xintercept = time_mid),linewidth=0.5,linetype = "dashed",color = "gray50")
  }
  if(label.sig & any.sig){
    probPlot <- probPlot +
      geom_label(data = sig.event,aes(x = time_mid,y = max.y,label = pvallab))
  }


  timeSeries <- plotSectionMeans(add.to.plot = ribbons,sig.event,time = timeMed,vals = valMed)+
    this_x_scale(name = x.axis.label, position = "top")+
    this_y_scale(name = y.axis.label, position = "right")+
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

#' Plot a histogram of null hypothesis and KDE results
#'
#' @param x a actR output object
#' @param pal Color brewer palette
#' @param h KDE smoothing factor, lower is smoother
#'
#' @return a ggplot object
#' @export
plotNullHistogram <- function(x,pal = "Paired",h = NA){
  nulls <- x$nullDetectionWithUncertainty[[1]]
  real <- x$eventDetectionWithUncertainty

  kd <- kdePval(nulls, real, h = h)

  hist <- ggplot()+
    actR_ggtheme()+
    geom_histogram(aes(x = nulls,y = after_stat(density),fill = "Null hypothesis results")) +
    geom_area(aes(x = kd$x,y = kd$y,fill = "KDE fit",color = "KDE fit"),alpha = 0.4) +
    geom_vline(aes(xintercept = real),color = "red") +
    geom_label(aes(x = real,y = max(kd$y) * 1.2),label = glue::glue("pval = {round(kd$pval,digits = 2)}")) +
    xlab("Fraction passed") +
    scale_fill_brewer(name = NULL,palette = pal) +
    scale_color_brewer(palette = pal) +
    theme(legend.position = c(.5,.8),legend.background = element_blank()) +
    guides(color = "none")

  return(hist)

}



#' Use a kde to estimate pvalue
#'
#' @param nulls a vector of null hypothesis results
#' @param real a single value for the real data
#' @param xmin minimum value for the distribution
#' @param xmax maximum value for the distribution
#' @param h smoothing factor
#'
#' @importFrom ks kde
#'
#' @return a list with pvalues and kde data
#' @export
kdePval <- function(nulls,real,h = NA,xmin = NA, xmax = NA){
  datarange <- range(c(nulls,real),na.rm = TRUE)
  span <- abs(diff(datarange))

  if(span == 0){
    if(real == 0){
      xmin <- 0
      xmax <- 1
      span <- 1
    }else{
      xmin <- real - 0.5
      xmax <- real + 0.5
      span <- 1
    }
  }

  if(all(is.na(h))){
    h <- .03 * span
  }

  if(all(is.na(xmin))){
    xmin <- min(datarange,na.rm = TRUE) - span*.05
  }

  if(all(is.na(xmax))){
    xmax <- max(datarange,na.rm = TRUE) + span*.05
  }

  #if xmin is less than 0, make symmetric
  if(xmin < 0){
    xmin <- -max(abs(c(xmin,xmax)))
    xmax <- max(abs(c(xmin,xmax)))
  }


#estimate kde
kd <- ks::kde(nulls,h = h * span,xmin = xmin,xmax = xmax,density = TRUE,gridsize = 1000)
if(real == xmax){
  real <- kd$eval.points[length(kd$eval.points)-1]
}

if(real == xmin){
  real <- kd$eval.points[2]
}

cdf <- cumsum(kd$estimate)/sum(kd$estimate)

pval <- 1-approx(kd$eval.points,y = cdf,xout = real)$y

out <- list(pval = pval,
            x = kd$eval.points,
            y = kd$estimate)

return(out)
}


#' Plot a timeseries result of sliding
#'
#' @param xlim limits on the plot
#' @param x sliding excursion output
#'
#' @return a plot object
#' @export
plotExcursionSliding <- function(x,xlim = c(12000,0)){

  # Make a basic plot
  p_values_all <- x$pvalue
  ages_to_test <- x$time_mid
  step <- abs(mean(diff(ages_to_test)))


# Set labels and title
xlabel <- paste('Age (',x$timeUnits[[1]],')',sep='')
ylabel <- paste(x$paleoData_proxy[[1]],' (',x$paleoData_units[[1]],')',sep='')
title  <- paste('Time series: ',x$archiveType[[1]],', ',x$dataSetName[[1]],', ',x$paleoData_TSid[[1]],
                ', lat=',x$geo_latitude[[1]],', lon=',x$geo_longitude[[1]],sep='')

# Top panel: plot the time series
par(mfrow=c(2,1))
plot(x$time[[1]],x$paleoData_values[[1]],type='b',xlim = xlim,xlab=xlabel,ylab=ylabel, main=title)

# Top panel: Use the p-values to shade the background
for (i in 1:length(p_values_all)) {
  if (!is.na(p_values_all[i])) {
    rect(xleft=ages_to_test[i]-(step/2),xright=ages_to_test[i]+(step/2),ybottom=par('usr')[3],ytop=par('usr')[4],
         col=adjustcolor('red',alpha=(1-p_values_all[i])/2),border='transparent')}
}
lines(x$time[[1]],x$paleoData_values[[1]],lwd=2) # Plot the record again, over top

# Bottom panel: plot the p-values
plot(ages_to_test,x$pvalue_positive,type='l',col = 'red',xlim=c(12000,0),ylim=c(1,0),xlab=xlabel,ylab="p-value",main='p-values (above in red, below in blue)')
lines(ages_to_test,x$pvalue_negative,col = 'blue')

abline(h=0.05,col='black',lty = 'dotted')

}



