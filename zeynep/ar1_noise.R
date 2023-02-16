#' @export
#' @name ar1.noise
#' @title Simulate n samples of an AR(1) process with given parameters
#' @description create synthetic AR(1) timeseries. Useful for null hypothesis testing
#' @param n number of samples (integer)
#' @param p number of variables (columns)
#' @param gamma the autocorrelation parameter of the model
#' @param std the standard deviation of the series
#'
#' @return a vector of synthetic values


ar1.noise <- function(n=100,ncol=1,gamma=0.7,sigma=1){
  ar1 = matrix(NA,nrow=n,ncol=ncol)
  for(j in 1:ncol){
    X <- arima.sim(model=list("ar"=gamma),n=n)
    ar1[,j] <- sigma*scale(X)
  }
  return(ar1)
}
