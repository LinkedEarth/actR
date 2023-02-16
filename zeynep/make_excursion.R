make_excursion <- function(length=1000, amp=1, del_amp=0, start=400, dur=200, shift=0, smooth=FALSE)
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
  return(series)
}
