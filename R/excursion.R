#' maximum number of consecutive TRUE values in a vector
#'
#' @param bool
#' @concept https://stackoverflow.com/questions/37345853/find-the-maximum-and-mean-length-of-the-consecutive-true-arguments
#' @return
#' @export
maxConsecutiveTRUE <- function(bool){

  consec <- rle(bool)

  # max consecutive values
  mv <- max(consec$lengths[consec$values==TRUE])
  if(!is.finite(mv)){
    mv <- 0
  }

  return(mv)
}


#' Title
#'
#' @param age
#' @param vals
#' @param event.yr
#' @param event.window
#' @param ref.window
#' @param n.consecutive
#' @param output.figure.path
#' @param sigNum
#'
#' @return
#' @export
#'
#' @examples
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
  if (preAVG + sigNum * preSD > postAVG + sigNum * postSD) {
    sd.hi = preSD
    avg.hi = preAVG
  } else {
    sd.hi = postSD
    avg.hi = postAVG
  }

  # Identify mean and sd used to test the low/negative anomaly threshold
  if (preAVG - sigNum * preSD < postAVG - sigNum * postSD) {
    sd.lo = preSD
    avg.lo = preAVG
  } else {
    sd.lo = postSD
    avg.lo = postAVG
  }

  # Identify points in event window that exceed the thresholds defined above

  aboveBool <- values[event.i] > avg.hi + sigNum * sd.hi
  belowBool <- values[event.i] < avg.lo - sigNum * sd.lo

  abovePts = which(aboveBool)
  belowPts = which(belowBool)

  # Determine whether there are any consecutive extreme points - this qualifies an event
   if (maxConsecutiveTRUE(aboveBool) >= n.consecutive | maxConsecutiveTRUE(belowBool) >= n.consecutive) {
    eventEX = TRUE

    # Check direction
    indsA = which(diff(abovePts) == 1)
    indsB = which(diff(belowPts) == 1)

    if (length(indsA) > 0 && length(indsB) > 0) {
      diffA = max(values[event.i[abovePts[c(indsA, indsA+1)]]] - (avg.hi + sigNum * sd.hi))
      diffB = max((avg.lo - sigNum * sd.lo) - values[event.i[belowPts[c(indsB, indsB+1)]]])
      dirEX = ifelse(diffA > diffB, 1, -1)
    } else {
      dirEX = ifelse(length(indsA) > 0, 1, -1)
    }

    # find event onset
    onsetEX = ifelse(dirEX == 1, age[event.i[abovePts[indsA[length(indsA)]+1]]], age[event.i[belowPts[indsB[length(indsB)]+1]]])

  }

  if (!is.na(output.figure.path)) {

  #plotExcursion() TBD
  }

  return(list(statusEX, eventEX, dirEX, onsetEX))

}
