#' Build a synthetic LiPD-timeseries-tibble with an excusrion
#'
#' @param length number of time steps in series
#' @param amplitude amplitude of excursion
#' @param delta.amplitude amplitude of change during excursion
#' @param start.time beginning time step of excursion
#' @param duration duration of excursion
#' @param time time vector over which to make the excursion
#' @param snr signal to noise ratio
#' @param sample.spacing.frac if < 1, fraction of the data to retain, randomly removing others
#' @param noise.ar1 AR1 coefficient for the noise
#' @param ... parameters passed to prepareInput
#'
#' @return ltt
#' @export
#'
makeExcursion <- function(time = NA,
                           length=1000,
                           amplitude=1,
                           delta.amplitude=0,
                           start.time=400,
                           duration=200,
                           snr = 1,
                           sample.spacing.frac = 1,
                           noise.ar1 = 0.5,
                           ...)
{

  if(all(is.na(time))){
  #Set the time axis
  time <- 1:length
  }

  length <- min(length,length(time))
  #Set the value axis
  values <- numeric(length)

  start <- which.min(abs(time - start.time))
  end <- which.min(abs(time - (start.time + duration)))

  #Compute the excursion line
  skew_axis <- seq(from=0,to=1, length.out=abs(end-start)+1)
  skew_values <- amplitude + delta.amplitude*skew_axis
  values[start:end] <- skew_values

  #Create series
  series <- values

  #subsample?
  if(sample.spacing.frac < 1){
    nsamps <- round(sample.spacing.frac * length)
    to.keep <- sort(sample(1:length,size = nsamps,replace = FALSE))
    series <- series[to.keep]
    time <- time[to.keep]
  }


  #Add noise
  if(is.finite(snr)){
    signal <- sd(series)
    noise <- simulateAutoCorrelatedUncertainty(sd = signal/snr,
                                               mean = 0,
                                               ar = noise.ar1,
                                               n = length(series))
    series <- noise + series
  }




  ltt <- prepareInput(vals = series,
                      time = time,
                      vals.variable.name = "Excursion test",
                      vals.units = "unitless",
                      dataset.name = "Synthetic excursion",
                      ...)

  return(ltt)
}

#' Build a synthetic LiPD-timeseries-tibble with a shift in mean
#'
#' @param length number of time steps
#' @param amp amplitude of shift
#' @param start starting time step of shift
#' @param dur duration of mean shift
#' @param shift value of time series before shift
#' @param ... values passed to prepareInput
#'
#' @return ltt
#' @export
#'
makeShift <- function(length, amp, start, dur=0, shift=0, ...)
{
  #Initialize our time axis
  time <- 0:(length-1)

  #Initialize values
  values <- integer(length)

  if (dur != 0)
  {
    jump_axis <- 0:dur
    jump_values <- (amp/dur)*jump_axis
    end <- dur + start
    values[start:end] <- jump_values
    values[end:length(values)] <- amp
  }
  else
  {
    values[start:length(values)] <- amp
  }

  series <- values + shift
  #time will be read as BP, so we will reverse the series
  ltt <- prepareInput(vals = rev(series), time = time,...)
  return(ltt)
}
