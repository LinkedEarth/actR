#' Summarize changepoint mean changes
#'
#' @param x Output of detectShiftCore()
#' @param alpha significance level of changepoints to consider
#' @import tibble
#'
#' @return a tibble of change point mean changes
#' @export
#'
#' @examples
summarizeChangepointMeanChanges <- function(x,alpha = 0.05){


  paramTib <- createTibbleFromParameterString(x$parameters[1])

  sig.event <- summarizeShiftSignificance(x$shiftDetection,
                                          alpha = alpha,
                                          paramTib = paramTib)

time <- apply(x$timeEns,1,median,na.rm = TRUE)
vals <- apply(x$valEns,1,median,na.rm = TRUE)

cpt <- sort(c(min(time,na.rm = TRUE),sig.event$time_start,max(time,na.rm = TRUE)))

sds <- ind <- means <- matrix(NA,nrow = length(time))

for(i in 1:(length(cpt)-1)){
  w <- which(dplyr::between(time,cpt[i],cpt[i+1]))
  ind[w] <- i
  means[w] <- mean(vals[w],na.rm = TRUE)
  sds[w] <- sd(vals[w],na.rm = TRUE)

}

tp <- tibble::tibble(time,
                     vals,
                     cpt_section = factor(ind),
                     cpt_mean = means) %>%
  dplyr::group_by(cpt_section) %>%
  dplyr::summarise(interval_time_min = min(time,na.rm = TRUE),
                   interval_time_max = max(time,na.rm = TRUE),
                   interval_mean = mean(cpt_mean,na.rm = TRUE))

tp$alpha = alpha

return(tp)
}


summarizeEventProbability <- function(exc.out,
                                      bin.step = 10,
                                      max.time = NA,
                                      min.time = NA){
  if(is.na(min.time)){
    min.time <- purrr::map_dbl(exc.out$time,min,na.rm=T) %>% median(na.rm = TRUE)
  }
  if(is.na(max.time)){
    max.time <- purrr::map_dbl(exc.out$time,max,na.rm=T) %>% median(na.rm = TRUE)
  }
  timeBins <- seq(min.time,max.time,by = bin.step)
  timeOut <- min.time+bin.step/2
  good.exc <- dplyr::filter(exc.out,eventDetected == TRUE) %>%
    dplyr::mutate(time_mid = (time_start + time_end)/2)


  eventSums <- eventsInWindow(good.exc$time_mid,start.vec = timeBins[-length(timeBins)],end.vec = timeBins[-1])

  out <- tibble::tibble(time_start = timeBins[-length(timeBins)],
                        time_end = timeBins[-1],
                        event_probability =  eventSums/exc.out$nEns[1]) %>%
          dplyr::mutate(time_mid = (time_start + time_end)/2)


  return(out)


}



eventsInWindow <- function(val,start.vec,end.vec){
  totalEvents <- purrr::map2(start.vec,end.vec,.f = ~ dplyr::between(val,.x,.y)) %>%
    unlist() %>%
    matrix(nrow = length(start.vec),byrow = TRUE) %>%
    rowSums()

  return(totalEvents)
}
