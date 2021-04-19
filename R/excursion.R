#' Maximum number of consecutive values in a vector
#'
#' @param vec a vector
#' @param val the value to search for consecutive values
#'
#' @concept https://stackoverflow.com/questions/37345853/find-the-maximum-and-mean-length-of-the-consecutive-true-arguments
#' @return the number of consecutive values. 0 means that value is not present in the vector
#' @export
#' @examples
#' maxConsecutive(c(1,2,7,5,7,7,7,3,4),val = 7)
maxConsecutive <- function(vec,val = TRUE){

  consec <- rle(vec)

  # max consecutive values
  mv <- max(consec$lengths[consec$values==val])
  if(!is.finite(mv)){
    mv <- 0
  }

  if(mv > 0){
    #get first indices for max
    fi <- min(which(consec$lengths == mv & consec$values == val))
    io <- seq(cumsum(consec$lengths)[fi]-(mv-1),cumsum(consec$lengths)[fi])
  }else{
    io <- NA
  }

  return(list(max = mv,index = io))
}


#' Detect an excursion in a timeseries
#'
#' @author Hannah Kolus
#' @author Nick McKay
#' @references Morrill
#' @param age a year or age vector
#' @param vals a vector of values the same length as age
#' @param event.yr the center of the proposed excursion window
#' @param event.window the width of the proposed excursion window
#' @param ref.window how many years to use as a reference before and after the event window
#' @param n.consecutive how many consecutive points are required for this to be considered an excursion? (default = 2)
#' @param output.figure.path path pointing to where should the output figure be saved? An NA will not produce a figure (default = NA)
#' @param sig.num how many standard deviations required outside the reference windows must be exceeded for this to be considered an excursion? (default = 2)
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
                           n.consecutive = 2,
                           output.figure.path = NA,
                           sig.num = 2) {
  ## Written by Hannah Kolus, 09/04/2018
  ## Determines whether an excursion event has occurred within the specified event window.
  ## Excursion events are defined as two consecutive values within the event window that
  ## are more extreme than the avg +/- X std of the reference windows.

  # yr.start:yr.end defines boundaries of analysis (i.e. both reference windows and the event window)
  yr.start = event.yr - event.window / 2 - ref.window
  yr.end = event.yr + event.window / 2 + ref.window

  # event.start:event.end defines the boundaries of the event
  event.start = event.yr - event.window / 2
  event.end = event.yr + event.window / 2

  analysis.i = which(age >= yr.start & age <= yr.end) # define analysis window indices

  age = age[analysis.i]
  vals = vals[analysis.i]

  # Detrend over analysis window
  a = predict(lm(vals ~ age))
  values = as.vector(vals - a)

  pre.i = which(age < event.start)                        # define pre-event (ref) window indices
  event.i = which(age >= event.start & age <= event.end)  # define event window indices
  post.i = which(age > event.end)                         # define post-event (ref) window indices

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
  # first look for positive excursions

  aboveBool <- values[event.i] > avg.hi + sig.num * sd.hi
  mctAbove <- maxConsecutive(aboveBool)
  abovePts = which(aboveBool)

  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (mctAbove$max >= n.consecutive) {
    valuesEX = age[mctAbove$index]

    outRowAbove <- tibble::tibble(excursionDetected = TRUE,
                             excursionDirection = 1,
                             excursionValues = list(valuesEX))
  }else{
    outRowAbove <- tibble::tibble(excursionDetected = FALSE,
                                  excursionDirection = 1,
                                  excursionValues = NA)
  }

  belowBool <- values[event.i] < avg.lo - sig.num * sd.lo
  mctBelow <- maxConsecutive(belowBool)

  # Determine whether there are any consecutive extreme points - this qualifies an event
  if (mctBelow$max >= n.consecutive) {
    # find event onset and termination
    valuesEX = age[mctAbove$index]

    outRowBelow <- tibble::tibble(excursionDetected = TRUE,
                                  excursionDirection = -1,
                                  excursionValues = list(valuesEX))
  }else{
    outRowBelow <- tibble::tibble(excursionDetected = FALSE,
                                  excursionDirection = -1,
                                  excursionValues = NA)
  }


  if (!is.na(output.figure.path)) {

    #plotExcursion() TBD
  }

  #prepare output tibblee
  out <- dplyr::bind_rows(outRowAbove,outRowBelow)

  return(out)

}
