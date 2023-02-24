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
#'
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


  #see if we got any results
  if(sum(!is.na(dataEst$eventDetected)) / nrow(dataEst) < 0.5){
    stop("Excursion detection couldn't be performed, probably because there were typically fewer than 'min.vals' poins than in the defined windows. Consider changing your parameters")
  }

  # now test null hypothesis
  nullHyp <- testNullHypothesis(time,
                                vals,
                                n.ens = n.ens,
                                surrogate.method = surrogate.method,
                                changeFun = detectExcursionCore,
                                how.long.prop = te,
                                mc.ens = null.hypothesis.n,
                                ...)

  nullEvents <- purrr::map_dbl(nullHyp,~ mean(.x$eventDetected,na.rm = TRUE)) %>%
    tibble::tibble(nulls = .)

  nullEcdf <- stats::ecdf(nullEvents$nulls)

  nullEventProb <- nullEvents %>%
    dplyr::summarize(qs = quantile(nulls,probs = null.quantiles))

  nullEventProb$clLevel <- paste0("cl",null.quantiles)

  nullLevels <- nullEventProb %>%
    tidyr::pivot_wider(values_from = qs,names_from = clLevel)

  eventSummary <- tibble::tibble(time_start = mean(dataEst$time_start,na.rm = TRUE),
                                 time_end = mean(dataEst$time_end,na.rm = TRUE),
                                 time_mid = mean(time_start,time_end),
                                 eventDetectionWithUncertainty = mean(dataEst$eventDetected,na.rm = TRUE),
                                 empirical_pvalue = 1-nullEcdf(eventDetectionWithUncertainty),
                                 eventDetection = list(dataEst),
                                 unc.prop.n = n.ens,
                                 null.hypothesis.n = null.hypothesis.n) %>%
    dplyr::bind_cols(nullLevels)

  paramTib <- createTibbleFromParameterString(as.character(dataEst$parameters[1]))

  eventSummary <- dplyr::bind_cols(eventSummary,paramTib,prepped)

# assign the appropriate class
  eventSummary <- new_excursion(eventSummary)

  return(eventSummary)

}


#' Detect excursion - core functionality
#'
#' @param time time vector of only the points in the window
#' @param vals value vector of only the points in the window
#' @param n.consecutive how many consecutive points are required for this to be considered an excursion? (default = 2)
#' @param exc.type Type of excursion to look for. "positive", "negative", "either" or "both" (default = "either")
#' @param min.vals Minimum number of values required in reference and event windows (default = 8)
#' @param na.rm Remove NAs? (default = TRUE)
#' @param sig.num how many standard deviations required outside the reference windows must be exceeded for this to be considered an excursion? (default = 2)
#' @param event.yr time at the center of the excursion window
#' @param event.window width (in time units) of the excursion window
#' @param ref.window width (in time units) of the reference windows
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

  # Detrend over analysis window
  a = predict(lm(vals ~ time))
  values = as.vector(vals - a)

  pre.i = which(time < event.start)                        # define pre-event (ref) window indices
  event.i = which(time >= event.start & time <= event.end)  # define event window indices
  post.i = which(time > event.end)                         # define post-event (ref) window indices

  #test for sufficient values in each window
  if(min(length(pre.i),length(event.i), length(post.i)) < min.vals){
    out <- tibble::tibble(time_start = event.start,
                          time_end = event.end,
                          eventDetected = NA,
                          eventProbability = NA,
                          time = list(time),
                          vals = list(values),
                          preMean = NA,
                          preSd = NA,
                          postMean = NA,
                          postSd = NA,
                          nExcursionVals = NA,
                          excursionMeanTime = NA,
                          excursionMaxSd = NA,
                          isExcursion = list(rep(NA,times = length(time))),
                          parameters = as.character(params))
    return(out)
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


  if(grepl(pattern = "either",x = exc.type,ignore.case = TRUE)){
    event <- any(c(aboveEvent,belowEvent))
    exc.ind <- c(exc.ind.above,exc.ind.below)
  }else if(grepl(pattern = "both",x = exc.type,ignore.case = TRUE)){
    event <- all(c(aboveEvent,belowEvent))
    exc.ind <- c(exc.ind.above,exc.ind.below)
  }else if(grepl(pattern = "pos",x = exc.type,ignore.case = TRUE)){
    event <- aboveEvent
    exc.ind <- exc.ind.above
  }else if(grepl(pattern = "neg",x = exc.type,ignore.case = TRUE)){
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



  out <- tibble::tibble(time_start = event.start,
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

  return(out)
}
