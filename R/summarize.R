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


  eventSums <- eventsInWindow(na.omit(good.exc$time_mid),start.vec = timeBins[-length(timeBins)],end.vec = timeBins[-1])

  out <- tibble::tibble(time_start = timeBins[-length(timeBins)],
                        time_end = timeBins[-1],
                        event_probability =  eventSums/exc.out$nEns[1]) %>%
          dplyr::mutate(time_mid = (time_start + time_end)/2)


  return(out)


}



eventsInWindow <- function(val,start.vec,end.vec){
  isna <- which(is.na(val))
  if(length(isna) > 0){
    if(length(isna) == length(val)){#if they're all NAs then 0 events
      return(matrix(0,nrow = length(start.vec)))
    }
    val <- val[-isna]#remove NAs
  }
  totalEvents <- purrr::map2(start.vec,end.vec,.f = ~ dplyr::between(val,.x,.y)) %>%
    unlist() %>%
    matrix(nrow = length(start.vec),byrow = TRUE) %>%
    rowSums()

  return(totalEvents)
}
