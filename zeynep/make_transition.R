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
  df = data.frame(time=time, values=series)
  return(df)
}
