summarizeEventProbability <- function(exc.out,
                                      bin.step = 10,
                                      max.age = NA,
                                      min.age = NA){
  if(is.na(min.age)){
    min.age <- purrr::map_dbl(exc.out$age,min,na.rm=T) %>% median(na.rm = TRUE)
  }
  if(is.na(max.age)){
    max.age <- purrr::map_dbl(exc.out$age,max,na.rm=T) %>% median(na.rm = TRUE)
  }
  ageBins <- seq(min.age,max.age,by = bin.step)
  ageOut <- min.age+bin.step/2
  good.exc <- dplyr::filter(exc.out,eventDetected == TRUE) %>%
    dplyr::mutate(time_mid = (time_start + time_end)/2)


  eventSums <- eventsInWindow(good.exc$time_mid,start.vec = ageBins[-length(ageBins)],end.vec = ageBins[-1])

  out <- tibble::tibble(time_start = ageBins[-length(ageBins)],
                        time_end = ageBins[-1],
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
