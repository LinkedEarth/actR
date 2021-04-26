detectExcursionSlidingWindow <- function(age,
                              vals,
                              event.window,
                              ref.window,
                              slide.step,
                              ...){

  #define the intervals over which to slide
  slide.min <- (min(age,na.rm = TRUE)+ref.window)+event.window/2
  slide.max <- (max(age,na.rm = TRUE)-ref.window)-event.window/2

  window.vec <- seq(slide.min,slide.max,by =slide.step)
  window.mid <- window.vec[-1]-event.window/2

  des <- function(ey,...){detectExcursion(event.yr = ey,...)}
  slid <- purrr::map_dfr(window.mid,
                         .f = des,
                         age = age,
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
#' @inheritParams propagateUncertainty
#' @param event.yr the center of the proposed excursion window
#' @param event.window the width of the proposed excursion window
#' @param ref.window how many years to use as a reference before and after the event window
#' @inheritDotParams propagateUncertainty
#' @param output.figure.path path pointing to where should the output figure be saved? An NA will not produce a figure (default = NA)
#'
#' @importFrom stats lm predict sd
#'
#' @return a tibble that describes the positive and negative excursion results
#' @export
detectExcursion = function(age,
                           vals,
                           event.yr,
                           event.window,
                           ref.window,
                           n.ens = 100,
                           output.figure.path = NA,
                           ...) {


  # yr.start:yr.end defines boundaries of analysis (i.e. both reference windows and the event window)
  yr.start = event.yr - event.window / 2 - ref.window
  yr.end = event.yr + event.window / 2 + ref.window

  # event.start:event.end defines the boundaries of the event
  event.start = event.yr - event.window / 2
  event.end = event.yr + event.window / 2

  analysis.i = which(age >= yr.start & age <= yr.end) # define analysis window indices

  age = age[analysis.i]
  vals = vals[analysis.i]

  # detect excursions while propagating age and data uncertainties
  dataEst <- propagateUncertainty(age,
                                  vals,
                                  n.ens = n.ens,
                                  changeFun = detectExcursionCore,
                                  event.start = event.start,
                                  event.end = event.end,
                                  ...)

  # now test null hypothesis
  nullHyp <- testNullHypothesis(age,
                                vals,
                                n.ens = n.ens,
                                method = "isopersistent",
                                changeFun = detectExcursionCore,
                                event.start = event.start,
                                event.end = event.end,
                                ...)

  eventSummary <- tibble::tibble(time_start = event.start,
                                 time_end = event.end,
                                 eventDetectionWithUncertainty = mean(dataEst$eventDetected),
                                 nullHypothesisDetection = mean(nullHyp$eventDetected),
                                 empiricalPvalue = nullHypothesisDetection/eventDetectionWithUncertainty,
                                 eventDetection = list(dataEst),
                                 nEns = n.ens)

  paramTib <- createTibbleFromParameterString(as.character(glue::glue("event.yr = {event.yr}, event.window = {event.window}, ref.window = {ref.window}, {dataEst$parameters}")))

  eventSummary <- dplyr::bind_cols(eventSummary,paramTib)

  eventSummary$empiricalPvalue[eventSummary$empiricalPvalue < 0] <- 0
  eventSummary$empiricalPvalue[eventSummary$empiricalPvalue > 1] <- 1


  if (!is.na(output.figure.path)) {
    #plotExcursion() TBD
  }

  return(eventSummary)

}


#' Detect excursion - core functionality
#'
#' @param age age vector of only the points in the window
#' @param vals value vector of only the points in the window
#' @param event.end the end of the event window
#' @param event.start the start of the event window
#' @param n.consecutive how many consecutive points are required for this to be considered an excursion? (default = 2)
#' @param exc.type Type of excursion to look for. "positive", "negativee", "either" or "both" (default = "either")
#' @param min.vals Minimum number of values required in reference and event windows (default = 8)
#' @param na.rm Remove NAs? (default = TRUE)
#' @param sig.num how many standard deviations required outside the reference windows must be exceeded for this to be considered an excursion? (default = 2)
#'
#' @return a tibble of results
#' @export
detectExcursionCore <- function(age,
                                vals,
                                event.start,
                                event.end,
                                sig.num = 2,
                                n.consecutive = 2,
                                exc.type = "either",
                                min.vals = 8,
                                na.rm = TRUE){

 #write parameters for export
  params = glue::glue("sig.num = {sig.num}, n.consecutive = {n.consecutive},exc.type = '{exc.type}', min.vals = {min.vals}, na.rm = {na.rm}")

  #removee NAs
  if(na.rm){
    good <- which(!is.na(age) & !is.na(vals))
    age <- age[good]
    vals <- vals[good]
  }

  # Detrend over analysis window
  a = predict(lm(vals ~ age))
  values = as.vector(vals - a)

  pre.i = which(age < event.start)                        # define pre-event (ref) window indices
  event.i = which(age >= event.start & age <= event.end)  # define event window indices
  post.i = which(age > event.end)                         # define post-event (ref) window indices

  #test for sufficient values in each window
  if(min(length(pre.i),length(event.i), length(post.i)) < min.vals){
    out <- tibble::tibble(time_start = event.start,
                          time_end = event.end,
                          eventDetected = NA,
                          eventProbability = NA,
                          age = list(age),
                          vals = list(values),
                          preMean = NA,
                          preSd = NA,
                          postMean = NA,
                          postSd = NA,
                          nExcusionVals = NA,
                          excursionMeanAge = NA,
                          excursionMaxSd = NA,
                          isExcursion = list(rep(NA,times = length(age))),
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
  exc.ind <- c() #setup excursion index

  # positive
  aboveBool <- values[event.i] > avg.hi + sig.num * sd.hi
  mctAbove <- maxConsecutive(aboveBool,gte = n.consecutive)

  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (mctAbove$max >= n.consecutive) {
    aboveEvent <- TRUE
    exc.ind <- c(exc.ind,mctAbove$index)
  }else{
    aboveEvent <- FALSE
  }

  belowBool <- values[event.i] < avg.lo - sig.num * sd.lo
  mctBelow <- maxConsecutive(belowBool,gte = n.consecutive)

  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (mctBelow$max >= n.consecutive) {
    belowEvent <- TRUE
    exc.ind <- c(exc.ind,mctBelow$index)
  }else{
    belowEvent <- FALSE
  }

  isExcursion <- seq_along(age) %in% event.i[exc.ind]

  if(any(isExcursion)){
    excursionMeanAge = mean(age[isExcursion])
    excursionMaxSd = max(abs(vals[isExcursion])/max(preSD,postSD))
  }else{
    excursionMeanAge = NA
    excursionMaxSd = NA
  }

  if(grepl(pattern = "either",x = exc.type,ignore.case = TRUE)){
    event <- any(c(aboveEvent,belowEvent))
  }else if(grepl(pattern = "both",x = exc.type,ignore.case = TRUE)){
    event <- all(c(aboveEvent,belowEvent))
  }else if(grepl(pattern = "pos",x = exc.type,ignore.case = TRUE)){
    event <- aboveEvent
  }else if(grepl(pattern = "neg",x = exc.type,ignore.case = TRUE)){
    event <- belowEvent
  }else{
    stop(glue::glue("exc.type = {exc.type} is not recognized"))
  }

  out <- tibble::tibble(time_start = event.start,
                        time_end = event.end,
                        eventDetected = event,
                        eventProbability = NA,
                        age = list(age),
                        vals = list(values),
                        preMean = preAVG,
                        preSd = preSD,
                        postMean = postAVG,
                        postSd = postSD,
                        nExcusionVals = sum(isExcursion),
                        excursionMeanAge = excursionMeanAge,
                        excursionMaxSd = excursionMaxSd,
                        isExcursion = list(isExcursion),
                        parameters = as.character(params))

  return(out)
}
