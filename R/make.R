#' Build a synthetic LiPD-timeseries-tibble with an excusrion
#'
#' @param length number of time steps in series
#' @param amp amplitude of excursion
#' @param del_amp amplitude of change during excursion
#' @param start begining time step of excursion
#' @param dur duration of excursion
#'
#' @return ltt
#' @export
#'
make_excursion <- function(length=1000, amp=1, del_amp=0, start=400, dur=200)
{
  #Set the time axis
  time <- 0:length

  #Set the value axis
  values <- numeric(length)

  #Compute the excursion line
  skew_axis <- seq(from=0,to=1, length.out=dur)
  skew_values <- amp + del_amp*skew_axis
  end <- start + (dur-1)
  values[start:end] <- skew_values

  #Create series
  series <- values

  #time will be read as BP, so we will reverse the series
  ltt <- prepareInput(vals = rev(series), time = 1:length(series))

  return(ltt)
}

#' Build a synthetic LiPD-timeseries-tibble with a shift in mean
#'
#' @param length number of time steps
#' @param amp amplitude of shift
#' @param start starting time step of shift
#' @param dur duration of mean shift
#' @param shift value of time series before shift
#' @param ...
#'
#' @return ltt
#' @export
#'
make_transition <- function(length, amp, start, dur=0, shift=0, ...)
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
  ltt <- prepareInput(vals = rev(series), time = time)
  return(ltt)
}
