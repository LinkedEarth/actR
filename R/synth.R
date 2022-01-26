# Function to create synthetic data for testing purposes

#' Linear ramp
#' @description linearRamp()  Create a linear ramp around the midpoint of a series. Slope can be adjusted as desired to create gradual or abrupt shifts.
#' @param lngth length of the time axis
#' @param width the width of the ramp. width  = 0 results in a step jump at the midpoint. width = lngth/2 (the largest alloweable value) results in a smooth linear ramp
#' @importFrom fBasics Ramp
#' @return values of the function (between 0 and 1)
#' @export
linearRamp <-function(lngth,width){
  if (width>lngth/2) {stop("Ramp width cannot exceed half of the length")}
  xs = seq(lngth)/lngth-1/2
  ramp = 4*fBasics::Ramp(xs,a = -width/lngth)
  ramp[xs>width/lngth] = 1.0
  return(ramp)
}

#' AR(1) noise
#' @description  samples from AR(1) process
#' @param lngth number of samples
#' @param sigma standard deviation of the Gaussian innovations
#' @param g lab-1 autocorrelation (short-term memory parameter)
#' @ImportFrom stats arima.sim
#' @return  lngth samples of the process (vector)
#' @export

ar1noise <- function(lngth,g,sigma){
  ar1=sigma/sqrt(1-g^2)*stats::arima.sim(model=list(g,0,0),n=lngth)
  return(ar1)
}

