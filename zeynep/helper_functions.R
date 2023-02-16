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

#noise function
ar1.noise <- function(n=100,ncol=1,gamma=0.7,sigma=1){
  ar1 = matrix(NA,nrow=n,ncol=ncol)
  for(j in 1:ncol){
    X <- arima.sim(model=list("ar"=gamma),n=n)
    ar1[,j] <- sigma*scale(X)
  }
  return(ar1)
}

make_excursion <- function(length=1000, amp=1, del_amp=0, start=400, dur=200, shift=0)
{
  #Set the time axis
  time <- 0:(length-1)
  
  #Set the value axis
  values <- numeric(length)
  
  #Compute the excursion line
  skew_axis <- seq(from=0,to=1, length.out=dur)
  skew_values <- amp + del_amp*skew_axis
  end <- start + (dur-1)
  values[start:end] <- skew_values
  
  #Create series
  series <- values
  series = data.frame(time=time, values=series)
  return(series)
}

