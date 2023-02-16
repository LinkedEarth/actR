#' Print Excursion
#'
#' @param x excursion object to print
#' @inheritDotParams summary.excursion
#'
#' @export
print.excursion <- function(x,...){
  summary.excursion(x,...)
}

#' Print shift
#'
#' @param x shift object to print
#' @inheritDotParams summary.shift
#'
#' @export
print.shift <- function(x,...){
  summary.shift(x,...)
}
