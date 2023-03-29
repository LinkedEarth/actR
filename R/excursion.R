detectExcursionSlidingWindow <- function(time,
                              vals,
                              event.window,
                              ref.window,
                              slide.step,
                              ...){

  #define the intervals over which to slide
  slide.min <- (min(time,na.rm = TRUE)+ref.window)+event.window/2
  slide.max <- (max(time,na.rm = TRUE)-ref.window)-event.window/2

  window.vec <- seq(slide.min,slide.max,by =slide.step)
  window.mid <- window.vec[-1]-event.window/2

  des <- function(ey,...){detectExcursion(event.yr = ey,...)}
  slid <- purrr::map_dfr(window.mid,
                         .f = des,
                         time = time,
                         vals = vals,
                         event.window = event.window,
                         ref.window = ref.window,
                         n.ens = 100,
                         ...)

  return(slid)

}

#helper function to allow pmapping over ts tibble rows
todfr <- function(...){
  l <- list(...)
  for(li in 1:length(l)){
    if(length(l[[li]]) > 1){
      l[[li]] <- list(l[[li]])
    }
    if(length(l[[li]]) == 0){
      l[[li]] <- NA
    }
  }
  return(tibble::as_tibble_row(l))
}

#' Detect an excursion in many timeseries
#'
#' @inheritParams detectExcursion
#' @inheritParams detectExcursionCore
#' @inheritParams propagateUncertainty
#' @inheritParams testNullHypothesis
#' @param seed Set a seed for reproducibility. By default it will use current time meaning it will not be reproducible.
#' @author Hannah Kolus
#' @author Nick McKay
#' @description Determines whether an excursion event has occurred within the specified event window for a lipd-ts-tibble of timeseries. Excursion events are defined as n.consecutive values within the event window that are more extreme than the avg +/- sig.num standard deviations of the reference windows.
#' @references Morrill
#'
#' @importFrom furrr future_pmap
#'
#' @return a tibble that describes the positive and negative excursion results
#' @export
detectMultipleExcursions <- function(ltt = NA,
                                     n.ens = 100,
                                     surrogate.method = "isospectral",
                                     null.hypothesis.n = 100,
                                     event.yr,
                                     event.window,
                                     ref.window,
                                     sig.num = 2,
                                     n.consecutive = 2,
                                     exc.type = "either",
                                     min.vals = 8,
                                     na.rm = TRUE,
                                     simulate.time.uncertainty = FALSE,
                                     simulate.paleo.uncertainty = FALSE,
                                     seed = as.integer(Sys.time())){


  #this requeires a lipd-tibble-ts with multiple timeseries
  if(all(is.na(ltt)) | !is.data.frame(ltt)){
    stop("detectMultipleExcursions requires lipd-ts-tibble input. See ?prepareInput for help")
  }

  if(nrow(ltt) < 2){
    stop("you must enter at least 2 rows in your lipd-ts-tibble to use detectMultipleExcursions")
  }

out <- furrr::future_pmap(ltt,\(...) detectExcursion(todfr(...),
                                                     n.ens = n.ens,
                                                     surrogate.method = surrogate.method,
                                                     null.hypothesis.n = null.hypothesis.n,
                                                     event.yr = event.yr,
                                                     event.window = event.window,
                                                     ref.window = ref.window,
                                                     sig.num = sig.num,
                                                     n.consecutive = n.consecutive,
                                                     exc.type = exc.type,
                                                     min.vals = min.vals,
                                                     na.rm = na.rm,
                                                     simulate.paleo.uncertainty = simulate.paleo.uncertainty,
                                                     simulate.time.uncertainty = simulate.time.uncertainty,
                                                     seed = seed,
                                                     progress = FALSE),
                          .progress = TRUE)


  return(out)

}


#' Detect an excursion in a timeseries
#'
#' @author Hannah Kolus
#' @author Nick McKay
#' @description Determines whether an excursion event has occurred within the specified event window. Excursion events are defined as n.consecutive values within the event window that are more extreme than the avg +/- sig.num standard deviations of the reference windows.
#' @references Morrill
#'
#' @inheritParams prepareInput
#' @inheritParams testNullHypothesis
#' @inheritParams detectShift
#' @inheritDotParams propagateUncertainty
#' @param output.figure.path path pointing to where should the output figure be saved? An NA will not produce a figure (default = NA)
#' @param pvalue.method method for estimating a pvalue. Options are "kde" (the default) which will use a KDE to estimate the pvalue relative to the null, or "ecdf" which will use an empirical cumulative distribution function.
#' @importFrom stats lm predict sd
#'
#' @return a tibble that describes the positive and negative excursion results
#' @export
detectExcursion = function(ltt = NA,
                           time = NA,
                           vals = NA,
                           time.variable.name = NA,
                           vals.variable.name = NA,
                           time.units = NA,
                           vals.units = NA,
                           dataset.name = NA,
                           n.ens = 100,
                           output.figure.path = NA,
                           surrogate.method = "isospectral",
                           null.hypothesis.n = 100,
                           null.quantiles = c(.95),
                           pvalue.method = "kde",
                           exc.type = "either",
                           ...) {

  #prep inputs
  prepped <- prepareInput(ltt = ltt,
                          time = time,
                          vals = vals,
                          time.variable.name = time.variable.name,
                          vals.variable.name = vals.variable.name,
                          time.units = time.units,
                          vals.units = vals.units,
                          dataset.name = dataset.name,
                          expecting.one.row = TRUE,
                          sort.by.time = TRUE,
                          remove.time.nas = TRUE)

  time <- prepped$time[[1]]
  vals <- prepped$paleoData_values[[1]]




  # detect excursions while propagating time and data uncertainties
  dataEst <- propagateUncertainty(time,
                                  vals,
                                  n.ens = n.ens,
                                  changeFun = detectExcursionCore,
                                  ...)

  # # cat("length(dataEst): ", length(dataEst), "\n")
  # # cat("names(dataEst): ", names(dataEst), "\n")
  # #
  # # cat("dim(dataEst$Orig): ", dim(dataEst$Orig), "\n")
  # # cat("names(dataEst$Orig): ", names(dataEst$Orig), "\n")
  # # cat("dataEst$Orig: ", dataEst$Orig$eventDetected, "\n")
  # #
  # # cat("dim(dataEst$All$pos): ", dim(dataEst$All$pos), "\n")
  # # cat("names(dataEst$All$pos): ", names(dataEst$All$pos), "\n")
  # # cat("dataEst$All$pos$eventDetected: ", dataEst$All$pos$eventDetected, "\n")
  #
  startTime <- mean(dataEst$either$time_start,na.rm = TRUE)
  endTime <- mean(dataEst$either$time_end,na.rm = TRUE)
  #
  # #see if we got any results
  #
  eventSummarySafe <- tibble::tibble(time_start = startTime,
                                 time_end = endTime,
                                 time_mid = mean(c(startTime,endTime)),
                                 eventDetectionWithUncertainty = NA,
                                 nullDetectionWithUncertainty = list(NA),
                                 empirical_pvalue = NA,
                                 eventDetection = dataEst[[eval(exc.type)]],
                                 unc.prop.n = n.ens,
                                 null.hypothesis.n = null.hypothesis.n)
  #
  eventSummarySafe[paste0("cl",null.quantiles)] <- NA
  #
  paramTib <- createTibbleFromParameterString(as.character(dataEst[[eval(exc.type)]]$parameters[1]))
  #
  eventSummarySafe <- dplyr::bind_cols(eventSummarySafe,paramTib,prepped)

  eventSummarySafe <- new_excursion(eventSummarySafe)
  if(sum(!is.na(dataEst[[eval(exc.type)]]$eventDetected)) / length(dataEst[[eval(exc.type)]]$eventDetected) < 0.5){#safely exit

    message("Returning safe summary, less than half of runs completed successfully")

  return(eventSummarySafe)

  }else{

    eventSummaryAll <- list()
  #
  #   # now test null hypothesis
    nullHyp <- testNullHypothesis(time,
                                  vals,
                                  n.ens = n.ens,
                                  surrogate.method = surrogate.method,
                                  changeFun = detectExcursionCore,
                                  mc.ens = null.hypothesis.n,
                                  ...)
  #
    exc.type.all <- c("either", "both", "pos", "neg")
  #
  #   dataEstIndices <- split(1:(n.ens*4), ceiling(seq_along(1:(n.ens*4)) / n.ens))
  #   nullIndices <- split(1:(n.ens*null.hypothesis.n), ceiling(seq_along(1:(n.ens*null.hypothesis.n)) / (n.ens)))
  #   countA <- 0
  #
  #
    for (zz in exc.type.all){

  #     countA <- countA + 1
  #     dataEstIndex <- dataEstIndices[[countA]]
  #
      cat(zz, "\n")
  #
  #     # cat("nullEvents: ", unlist(lapply(nullHyp, function(x) x$Orig$eventDetected[dataEstIndex])), "\n")
  #     # cat("numnullevents: ", length(unlist(lapply(nullHyp, function(x) x$Orig$eventDetected[dataEstIndex]))), "\n")
  #     #
  #     # allNULLS <- unlist(lapply(nullHyp, function(x) x$Orig$eventDetected[dataEstIndex]))
  #     #
  #     # nullEvents <- list()
  #     # for (i in 1:length(nullIndices)){
  #     #   iii <- nullIndices[[i]]
  #     #   nullEvents[[i]] <- mean(iii)
  #     # }
  #     # nullEvents <- unlist(nullEvents)

      cat(nullHyp[[1]][[eval(zz)]]$eventDetected, "\n")
  #
      nullEvents <- purrr::map_dbl(nullHyp,~ mean(.x[[eval(zz)]]$eventDetected,na.rm = TRUE)) %>%
        tibble::tibble(nulls = .)
      if(all(!is.finite(nullEvents$nulls))){
        eventSummaryAll[[eval(zz)]] <- eventSummarySafe
        next
      }
      cat(nullEvents$nulls, "\n")
  #
      eventDetectionWithUncertainty <-  mean(dataEst[[eval(zz)]]$eventDetected,na.rm = TRUE)
  #
  #     #cat(dataEst$Orig$eventDetected[dataEstIndex], "\n")
  #
  #     #cat("nullevents: ", nullEvents$nulls, "\n")
  #
  #
      if(pvalue.method == "ecdf"){

        nullEcdf <- stats::ecdf(nullEvents$nulls)
        pval <- 1-nullEcdf(eventDetectionWithUncertainty)

      }else if(pvalue.method == "kde"){
        pval <- kdePval(nullEvents$nulls,eventDetectionWithUncertainty)$pval
      }else{
        stop("pvalue.method must be 'ecdf' or 'kde'")
      }
  #
      nullEventProb <- nullEvents %>%
        dplyr::summarize(qs = quantile(nulls,probs = null.quantiles,na.rm = TRUE))

      nullEventProb$clLevel <- paste0("cl",null.quantiles)

      nullLevels <- nullEventProb %>%
        tidyr::pivot_wider(values_from = qs,names_from = clLevel)

      eventSummary <- tibble::tibble(time_start = startTime,
                                     time_end = endTime,
                                     time_mid = mean(c(startTime,endTime)),
                                     eventDetectionWithUncertainty = eventDetectionWithUncertainty,
                                     nullDetectionWithUncertainty = nullEvents$nulls,
                                     empirical_pvalue = pval,
                                     eventDetection = dataEst[[eval(zz)]],
                                     unc.prop.n = n.ens,
                                     null.hypothesis.n = null.hypothesis.n) %>%
        dplyr::bind_cols(nullLevels)

      paramTib <- createTibbleFromParameterString(as.character(dataEst[[eval(zz)]]$parameters[1]))

      eventSummary <- dplyr::bind_cols(eventSummary,paramTib,prepped)
  #
  #     # assign the appropriate class
      eventSummary <- new_excursion(eventSummary)
  #
      eventSummaryAll[[eval(zz)]]<-eventSummary
  #
    }
  #
  #
    # if (length(eventSummaryAll)<1){
    #   returns <- NULL
    # }else{
    #   returns <- eventSummaryAll
    # }
  #
  }
  #
  # returns <- list(returns, nullEvents)

    # returns <- list(dataEst, nullHyp)

  return(eventSummaryAll)

}


#' Detect excursion - core functionality
#'
#' @param time time vector of only the points in the window
#' @param vals value vector of only the points in the window
#' @param n.consecutive how many consecutive points are required for this to be considered an excursion? (default = 2)
#' @param exc.type Type of excursion to look for. "positive", "negative", "either" or "both" (default = "either")
#' @param min.vals Minimum effective sample size (adjusted by autocorrelation) required in reference and event windows (default = 4)
#' @param na.rm Remove NAs? (default = TRUE)
#' @param sig.num how many standard deviations required outside the reference windows must be exceeded for this to be considered an excursion? (default = 2)
#' @param event.yr time at the center of the excursion window
#' @param event.window width (in time units) of the excursion window
#' @param ref.window width (in time units) of the reference windows
#' @param adjust.n.for.autocorrelation Adjust the min.vals check to account for autocorrelation?
#'
#' @return a tibble of results
#' @export
detectExcursionCore <- function(time,
                                vals,
                                event.yr,
                                event.window,
                                ref.window,
                                sig.num = 2,
                                n.consecutive = 2,
                                exc.type = "either",
                                min.vals = 8,
                                adjust.n.for.autocorrelation = FALSE,
                                na.rm = TRUE){

 #write parameters for export
  params = glue::glue("event.yr = {event.yr}, event.window = {event.window}, ref.window = {ref.window}, sig.num = {sig.num}, n.consecutive = {n.consecutive},exc.type = '{exc.type}', min.vals = {min.vals}, na.rm = {na.rm}")

  #removee NAs
  if(na.rm){
    good <- which(!is.na(time) & !is.na(vals))
    time <- time[good]
    vals <- vals[good]
  }

  # yr.start:yr.end defines boundaries of analysis (i.e. both reference windows and the event window)
  yr.start = event.yr - event.window / 2 - ref.window
  yr.end = event.yr + event.window / 2 + ref.window

  # event.start:event.end defines the boundaries of the event
  event.start = event.yr - event.window / 2
  event.end = event.yr + event.window / 2

  analysis.i = which(time >= yr.start & time <= yr.end) # define analysis window indices

  time = time[analysis.i]
  vals = vals[analysis.i]

  #create safe output in case of error
  safeOut <- tibble::tibble(time_start = event.start,
                        time_end = event.end,
                        eventDetected = NA,
                        eventProbability = NA,
                        time = list(time),
                        vals = list(vals),
                        preMean = NA,
                        preSd = NA,
                        postMean = NA,
                        postSd = NA,
                        nExcursionVals = NA,
                        excursionMeanTime = NA,
                        excursionMaxSd = NA,
                        isExcursion = list(rep(NA,times = length(time))),
                        parameters = as.character(params))

  if(length(vals) <= 3){
    return(safeOut)
  }
  # Detrend over analysis window
  a <- try(predict(lm(vals ~ time)),silent = TRUE)
  if(is(a,"try-error")){
    return(safeOut)
  }
  values = as.vector(vals - a)

  pre.i = which(time < event.start)                        # define pre-event (ref) window indices
  event.i = which(time >= event.start & time <= event.end)  # define event window indices
  post.i = which(time > event.end)                         # define post-event (ref) window indices


  if(adjust.n.for.autocorrelation){
    ar <- arCumulative(vals)
    effective.n.adjustment <- (1)/(1+ar) #https://andrewcharlesjones.github.io/journal/21-effective-sample-size.html
  }else{
    effective.n.adjustment <- 1
  }

  #test for sufficient values in each window
  if(min( length(pre.i)*effective.n.adjustment, length(event.i), length(post.i)*effective.n.adjustment )  < min.vals){
    warning("insufficient minimum values")
    return(safeOut)
  }


  # Calculate the avg and sd for pre-event ref window, excluding the most extreme data point
  preAVG = mean(values[pre.i])
  extremeInd = which(max(abs(preAVG - values[pre.i])) == abs(preAVG - values[pre.i]))
  preAVG = mean(values[pre.i[-extremeInd[1]]])
  preSD = sd(values[pre.i[-extremeInd[1]]])

  # Calculate the avg and sd for post-event ref window, excluding the most extreme data point
  postAVG = mean(values[post.i])
  extremeInd = which(max(abs(postAVG - values[post.i])) == abs(postAVG - values[post.i]))
  postAVG = mean(values[post.i[-extremeInd[1]]])
  postSD = sd(values[post.i[-extremeInd[1]]])

  # Identify mean and sd used to test the high/positive anomaly threshold
  if (preAVG + sig.num * preSD > postAVG + sig.num * postSD) {
    sd.hi = preSD
    avg.hi = preAVG
  } else {
    sd.hi = postSD
    avg.hi = postAVG
  }

  # Identify mean and sd used to test the low/negative anomaly threshold
  if (preAVG - sig.num * preSD < postAVG - sig.num * postSD) {
    sd.lo = preSD
    avg.lo = preAVG
  } else {
    sd.lo = postSD
    avg.lo = postAVG
  }



  # Identify points in event window that exceed the thresholds defined above
  exc.ind.above <- c() #setup excursion index
  exc.ind.below <- c() #setup excursion index

  # positive
  aboveBool <- values[event.i] > avg.hi + sig.num * sd.hi
  mctAbove <- maxConsecutive(aboveBool,gte = n.consecutive)

  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (mctAbove$max >= n.consecutive) {
    aboveEvent <- TRUE
    exc.ind.above <- c(exc.ind.above,mctAbove$index)
  }else{
    aboveEvent <- FALSE
  }

  belowBool <- values[event.i] < avg.lo - sig.num * sd.lo
  mctBelow <- maxConsecutive(belowBool,gte = n.consecutive)

  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (mctBelow$max >= n.consecutive) {
    belowEvent <- TRUE
    exc.ind.below <- c(exc.ind.below,mctBelow$index)
  }else{
    belowEvent <- FALSE
  }

  out <- list()


  exc.type.all <- c("either", "both", "pos", "neg")


  for (zz in exc.type.all){
    if(grepl(pattern = "either",x = zz,ignore.case = TRUE)){
      event <- any(c(aboveEvent,belowEvent))
      exc.ind <- c(exc.ind.above,exc.ind.below)
    }else if(grepl(pattern = "both",x = zz,ignore.case = TRUE)){
      event <- all(c(aboveEvent,belowEvent))
      exc.ind <- c(exc.ind.above,exc.ind.below)
    }else if(grepl(pattern = "pos",x = zz,ignore.case = TRUE)){
      event <- aboveEvent
      exc.ind <- exc.ind.above
    }else if(grepl(pattern = "neg",x = zz,ignore.case = TRUE)){
      event <- belowEvent
      exc.ind <- exc.ind.below
    }else{
      stop(glue::glue("exc.type = {exc.type} is not recognized"))
    }

    isExcursion <- seq_along(time) %in% event.i[exc.ind]

    if(any(isExcursion)){
      excursionMeanTime = mean(time[isExcursion])
      excursionMaxSd = max(abs(vals[isExcursion])/max(preSD,postSD))
    }else{
      excursionMeanTime = NA
      excursionMaxSd = NA
    }

    out[[eval(zz)]] <- tibble::tibble(time_start = event.start,
                          time_end = event.end,
                          eventDetected = event,
                          eventProbability = NA,
                          time = list(time),
                          vals = list(values),
                          preMean = preAVG,
                          preSd = preSD,
                          postMean = postAVG,
                          postSd = postSD,
                          nExcursionVals = sum(isExcursion),
                          excursionMeanTime = excursionMeanTime,
                          excursionMaxSd = excursionMaxSd,
                          isExcursion = list(isExcursion),
                          parameters = as.character(params)) %>%
      new_excursionCore()
  }

  #returns<-list("Orig" = out[[eval(exc.type)]], "All" = out)

  return(out)
}

arCumulative <- function(x){
  a <- acf(x,plot = FALSE)
  sig <- qnorm((1 + .95)/2)/sqrt(a$n.used)
  wa <- c()
  ari <- 2
  while(a$acf[ari] > sig){
    wa <- c(wa,ari)
    ari <- ari + 1
    if(ari > length(a$acf)){break}
  }

  if(length(wa) == 0){
    arc <- 0
  }else{
    arc <- sum(a$acf[wa])
  }

  return(arc)
}
