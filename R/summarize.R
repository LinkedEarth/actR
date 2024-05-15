#' Summarize changepoint mean changes
#'
#' @param x Output of detectShiftCore
#' @param alpha significance level of changepoints to consider
#' @import tibble
#'
#' @return a tibble of change point mean changes
#' @export
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


#' Summarize the shift detection results using a uniform bin step
#'
#' @param exc.out   shiftCore (tbl_df) class result of propagateUncertainty()
#' @param bin.vec   Vector by which to bin the data
#' @param bin.step  Number of years to group shifts by
#' @param max.time  Maximum age of bin.vec to summarize
#' @param min.time  Minimum age of bin.vec to summarize
#' @param shift.type Type of excursion to look for. "positive", "negative", "either" or "both" (default = "either")
#'
#' @return a tibble of event probability for each bin step
#' @export
summarizeEventProbability <- function(exc.out,
                                      bin.vec = NA,
                                      bin.step = 10,
                                      max.time = NA,
                                      min.time = NA,
                                      shift.type = "either",
                                      time.dir = "retrograde"
                                      ){
  if(is.na(min.time)){
    min.time <- purrr::map_dbl(exc.out$time,min,na.rm=T) %>% median(na.rm = TRUE)
  }
  if(is.na(max.time)){
    max.time <- purrr::map_dbl(exc.out$time,max,na.rm=T) %>% median(na.rm = TRUE)
  }

  if(!all(is.na(bin.vec))){#make one using bin.step
    timeBins <- bin.vec
    bin.step <- median(abs(diff(timeBins)),na.rm = TRUE)
    min.time <- min(timeBins)
    max.time <- max(timeBins)
  }else{
    timeBins <- seq(min.time,max.time,by = bin.step)
  }

  #timeOut <- min.time+bin.step/2
  good.exc <- dplyr::filter(exc.out,eventDetected == TRUE) %>%
    dplyr::mutate(time_mid = (time_start + time_end)/2)


  eventSumsEither   <- eventsInWindow(na.omit(good.exc$time_mid),
                                      start.vec = timeBins[-length(timeBins)],end.vec = timeBins[-1])

  is.cpt.var <- grepl("cpt.var",exc.out$parameters[[1]])


  if(time.dir == "retrograde" & is.cpt.var){
    filtPos <- dplyr::filter(good.exc,delta_sd < 0)
  }else if(time.dir == "prograde" & is.cpt.var){
    filtPos <- dplyr::filter(good.exc,delta_sd >= 0)
  }else if(time.dir == "retrograde" & !is.cpt.var){
    filtPos <- dplyr::filter(good.exc,delta_mean < 0)
  }else if(time.dir == "prograde" & !is.cpt.var){
    filtPos <-dplyr::filter(good.exc,delta_mean >= 0)
  }

  if(time.dir == "retrograde" & is.cpt.var){
    filtNeg <- dplyr::filter(good.exc,delta_sd >= 0)
  }else if(time.dir == "prograde" & is.cpt.var){
    filtNeg <- dplyr::filter(good.exc,delta_sd < 0)
  }else if(time.dir == "retrograde" & !is.cpt.var){
    filtNeg <- dplyr::filter(good.exc,delta_mean >= 0)
  }else if(time.dir == "prograde" & !is.cpt.var){
    filtNeg <-dplyr::filter(good.exc,delta_mean < 0)
  }


  eventSumsPositive <- eventsInWindow(na.omit(filtPos$time_mid),start.vec = timeBins[-length(timeBins)],end.vec = timeBins[-1])
  eventSumsNegative <- eventsInWindow(na.omit(filtNeg$time_mid),start.vec = timeBins[-length(timeBins)],end.vec = timeBins[-1])
  eventSumsBoth     <- apply(matrix(c(eventSumsPositive,eventSumsNegative),length(eventSumsNegative)),1,min)

  #now pick the one that was chosen

  if(grepl(pattern = "either",x = shift.type,ignore.case = TRUE)){
    event <- eventSumsEither
  }else if(grepl(pattern = "both",x = shift.type,ignore.case = TRUE)){
    event <- eventSumsBoth
  }else if(grepl(pattern = "pos",x = shift.type,ignore.case = TRUE)){
    event <- eventSumsPositive
  }else if(grepl(pattern = "neg",x = shift.type,ignore.case = TRUE)){
    event <- eventSumsNegative
  }else{
    stop(glue::glue("shift.type = {shift.type} is not recognized"))
  }

  #Create output tibble
  out <- tibble::tibble(time_start = timeBins[-length(timeBins)],
                        time_end = timeBins[-1]
                        ) %>%
          dplyr::mutate(time_mid = (time_start + time_end)/2)

  out$event_probability          =  event/exc.out$nEns[1]
  out$event_probability_either   =  eventSumsEither/exc.out$nEns[1]
  out$event_probability_positive =  eventSumsPositive/exc.out$nEns[1]
  out$event_probability_negative =  eventSumsNegative/exc.out$nEns[1]
  out$event_probability_both     =  eventSumsBoth/exc.out$nEns[1]
  out$event_probability_both     =  eventSumsBoth/exc.out$nEns[1]

  #add delta for each shift
  out$deltas  <- NA

  deltas <- good.exc %>%
         mutate(bin = cut(time_mid, breaks = timeBins, labels = FALSE)) %>%
         group_by(bin)

  if(grepl("cpt.var",exc.out$parameters[[1]])){
    deltas <- deltas %>% summarise(mean_value = median(delta_sd))
  }else{
    deltas <- deltas %>% summarise(mean_value = median(delta_mean))
  }
  deltas <- deltas %>% filter(!is.na(bin))

  out$deltas[deltas$bin] <- deltas$mean_value*-1


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


rounder <- function(x){
  round(x,max(0,ceiling(1-log10(x))))
}


summarizeParams <- function(x){
  if(length(unique(x)) == 1){
    return(as.character(unique(x)))
  }else{
    if(all(is.numeric(x))){
      if(length(unique(x)) > 10){
        return(glue::glue("{rounder(mean(x))} ± {rounder(sd(x))}"))
      }else{
        return(paste(sort(unique(x)),collapse = ", "))
      }
    }else{
      return(paste(unique(x),collapse = ", "))
    }
  }
}
