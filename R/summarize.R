summarizeEventProbability <- function(exc.out,bin.step = 1){
  minAge <- purrr::map_dbl(exc.out$age,min,na.rm=T) %>% min(na.rm = TRUE)
  maxAge <- purrr::map_dbl(exc.out$age,max,na.rm=T) %>% max(na.rm = TRUE)
  ageBins <- seq(minAge,maxAge,by = bin.step)
  ageOut <- minAge+bin.step/2
  good.exc <- dplyr::filter(exc.out,eventDetected == TRUE) %>%
    dplyr::mutate(time_mid = (time_start + time_end)/2)


  eventSums <- eventsInWindow(good.exc$time_mid,start.vec = ageBins[-length(ageBins)],end.vec = ageBins[-1])

  out <- tibble::tibble(bin_start = ageBins[-length(ageBins)],
                        bin_end = ageBins[-1],
                        event_probability =  eventSums/exc.out$nEns[1]) %>%
          dplyr::mutate(bin_mid = (bin_start + bin_end)/2)


  return(out)


}



eventsInWindow <- function(val,start.vec,end.vec){
  totalEvents <- purrr::map2(start.vec,end.vec,.f = ~ dplyr::between(val,.x,.y)) %>%
    unlist() %>%
    matrix(nrow = length(start.vec),byrow = TRUE) %>%
    rowSums()

  return(totalEvents)
}
