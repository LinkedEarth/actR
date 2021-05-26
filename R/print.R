#' Print Excursion
#'
#' @param x excursion object to print
#' @inheritDotParams print
#'
#' @export
print.excursion <- function(x,...){
  summary.excursion(x,...)
}

#' Print shift
#'
#' @param x shift object to print
#' @inheritDotParams print
#'
#' @export
print.shift <- function(x,...){
  summary.shift(x,...)
}
