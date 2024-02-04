
#' Use the write function from rEDM depending on the version
#'
#' @param ... Arguments to pass to surrogate
#'
#' @return Surrogate output
#' @export
surrogateDataFun <- function(...){
  if(packageVersion("rEDM") < "1.15"){
    out <- rEDM::make_surrogate_data(...)
  }else{
    out <- rEDM::SurrogateData(...)
  }
  return(out)
}


#' Explore uncertainty space on abscissa or ordinate and propagate of any of actR's change functions
#'
#' @param time a time vector, or matrix of time ensemble members (ensembles in columns)
#' @param vals a values vector, or matrix of values ensemble members (ensembles in columns)
#' @param changeFun the change function to across which to propagate
#' @param simulate.time.uncertainty TRUE or FALSE. If an ensemble is not included, do you want to simulate time ensembles (default = TRUE)
#' @param simulate.paleo.uncertainty TRUE or FALSE. If an ensemble is not included, do you want to simulate paleo ensembles (default = TRUE)
#' @param n.ens How many ensembles to use for error propagation? (default = 100)
#' @param bam.model BAM Model parameters to use if simulating time uncertainty (default = list(ns = n.ens, name = "bernoulli", param = 0.05), paleo.uncertainty = sd(vals,na.rm = TRUE)))
#' @param paleo.uncertainty Uncertainty to use if modelling uncertainty for paleo values. (default = sd(vals,na.rm = TRUE)/2)
#' @param paleo.ar1 Autocorrelation coefficient to use for modelling uncertainty on paleoData, what fraction of the uncertainties are autocorrelated? (default = sqrt(0.5); or 50% autocorrelated uncertainty)
#' @param paleo.arima.order Order to use for ARIMA model used in modelling uncertainty on paleoDat (default = c(1,0,0))
#' @param summarize Boolean. Summarize the output? Or return all the ensembles?
#' @param ... arguments to pass to pass to changeFun
#' @param seed set a seed for reproducibility
#' @param progress show null hypothesis testing progress bar?
#'
#' @return a propagated uncertainty tibble
#' @export
propagateUncertainty <- function(time,
                                 vals,
                                 changeFun,
                                 simulate.time.uncertainty = TRUE,
                                 simulate.paleo.uncertainty = TRUE,
                                 n.ens = 100,
                                 bam.model = list(ns = n.ens, name = "bernoulli", param = 0.05),
                                 paleo.uncertainty = sd(vals,na.rm = TRUE)/2,
                                 paleo.ar1 = sqrt(0.5),
                                 paleo.arima.order = c(1,0,0),
                                 summarize = FALSE,
                                 seed = round(sum(time,na.rm = TRUE)),
                                 progress = TRUE,
                                 ...){




  #check inputs
  if(max(c(NCOL(time),NCOL(vals))) > 1 & n.ens <= 1){#at least one is an ensemble
    stop("To simulate uncertainty with an ensemble, increasee n.ens to more than 1 (probably more than 50 at the minimum)")
  }

  nca <- NCOL(time)

  #set a seed.
  if(any(is.na(seed))){
    seed <- sample.int(1000,1)
  }
  try(set.seed(seed),silent = TRUE)

  #Prepare time ensemble
  if(nca == 1){#then it's not an ensemble
    #create ensemble?
    if(simulate.time.uncertainty){
      timeMat <- geoChronR::simulateBam(X = matrix(1,nrow = length(time)),
                                       t = as.matrix(time),
                                       model = bam.model,
                                       ageEnsOut = TRUE)$ageEns
    }else{#replicate times up to n.ens
      timeMat <- matrix(rep(time,n.ens),ncol = n.ens,nrow = NROW(time))
    }
  }else{
    if(nca >= n.ens){
      timeMat <- time[,sample(seq_len(nca),size = n.ens,replace = FALSE)]
    }else if(nca < n.ens){
      timeMat <- time[,sample(seq_len(nca),size = n.ens,replace = TRUE)]
    }
  }

  #make into a list for purrr
  timeList <- purrr::array_branch(timeMat,margin = 2)


  #Now prep paleodata
  ncp <- NCOL(vals)

  if(ncp == 1){#then it's not an ensemble
    #create ensemble?
    if(simulate.paleo.uncertainty){
      paleoList <- purrr::rerun(vals + simulateAutoCorrelatedUncertainty(sd = paleo.uncertainty,
                                                                         n = NROW(vals),
                                                                         ar = paleo.ar1,
                                                                         arima.order = paleo.arima.order),.n = n.ens)


    }else{#replicate times up to n.ens
      paleoList <- matrix(rep(vals,n.ens),ncol = n.ens,nrow = NROW(vals)) %>%
        purrr::array_branch(margin = 2)

    }
  }else{
    if(ncp >= n.ens){
      paleoMat <- vals[,sample(seq_len(ncp),size = n.ens,replace = FALSE)]
    }else if(ncp < n.ens){
      paleoMat <- vals[,sample(seq_len(ncp),size = n.ens,replace = TRUE)]
    }
    paleoList <- purrr::array_branch(paleoMat,margin = 2)
  }

  # check to see if ... includes vectors of parameters
  dots <- list(...)

  dl <- purrr::map_dbl(dots,length)
  if(any(dl > 1)){#then we need to sample over
    #turn params into vectors that are n.ens long...
    dotsLong <- vector(mode = "list",length = length(dots))
    dotsSummaryString <- c()
    for(d in 1:length(dots)){
      if(length(dots[[d]]) == 1){
        dotsLong[[d]] <- rep(dots[[d]],n.ens)
        dotsSummaryString[d] <- glue::glue("{dots[[d]]}")
      }else if(length(dots[[d]]) <= n.ens){
        dotsLong[[d]] <- sample(dots[[d]],size = n.ens,replace = FALSE)
        dotsSummaryString[d] <- glue::glue("{round(mean(dotsLong[[d]]))} ± {round(sd(dotsLong[[d]]))}")
      }else{
        dotsLong[[d]] <- sample(dots[[d]],size = n.ens,replace = TRUE)
        dotsSummaryString[d] <- glue::glue("{round(mean(dotsLong[[d]]))} ± {round(sd(dotsLong[[d]]))}")
      }
    }
    names(dotsLong) <- names(dots)
    tomaplist <- append(list(time = timeList,vals = paleoList),dotsLong)

    propagated <- purrr::pmap_dfr(tomaplist,changeFun)
  }else{
    propagated <- purrr::map2_dfr(timeList,paleoList,changeFun,...)
  }

  propagated$nEns <- n.ens

  #not sure I want this
  if(summarize){
    propagated <- propagated %>%
      dplyr::group_by(time_start,time_end) %>%
      dplyr::summarise(meanDetected = mean(eventDetected),
                       meanProbability = mean(eventProbability),
                       nEns = n.ens,
                       parameters = unique(parameters))
  }

  return(propagated)



}


#' Detect excursions in synthetic datasets that mimic a real on
#'
#' @inheritParams propagateUncertainty
#' @param mc.ens How many Monte Carlo simulations to use for null hypothesis testing
#' @param surrogate.method What method to use to generage surrogate data for hypothesis testing? Options include: \itemize{
#' \item 'isospectral': (Default) Following Ebisuzaki (1997), generate surrogates by scrambling the phases of the data while preserving their power spectrum. This uses the To generate these “isospectral” surrogates. Uses the rEDM::make_surrogate_data() or rEDM::SurrogateData() function depending on version
#' \item 'isopersistent':  Generates surrogates by simulating from an autoregressive process of order 1 (AR(1)), which has been fit to the data. Uses the geoChronR::createSyntheticTimeseries() function
#' \item 'shuffle': Randomly shuffles the data to create surrogates. Uses the rEDM::make_surrogate_data() or rEDM::SurrogateData() function depending on version
#' }
#' @inheritDotParams propagateUncertainty
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @return a tibble that reports the positivity rate in the synthetics
#' @export
testNullHypothesis <- function(time,
                               vals,
                               changeFun,
                               n.ens = 100,
                               mc.ens = 100,
                               surrogate.method = "isospectral",
                               seed = round(sum(vals,na.rm=FALSE)),
                               progress = TRUE,
                               ...) {


  #prep values for surrogates

  #set a seed.
  if(all(is.na(seed))){
    seed <- sample.int(1000,1)
  }
  set.seed(seed)


  ncp <- NCOL(vals)

  if(ncp == 1){
    vm <- matrix(rep(vals,times = mc.ens),ncol = mc.ens,byrow = FALSE)
    valList <- purrr::array_branch(vm,margin = 2)
  }else{
    valList <- purrr::array_branch(vals,margin = 2)
    if(length(valList) >= n.ens){
      valList <- valList[sample(seq_along(valList),size = mc.ens,replace = FALSE)]
    }else{
      valList <- valList[sample(seq_along(valList),size = mc.ens,replace = TRUE)]
    }
  }


  # generate some surrogates

  if(grepl(surrogate.method,pattern = "persis",ignore.case = T)){
    #surVals <- geoChronR::ar1Surrogates(time = time,vals = vals,detrend = TRUE,method = "redfit",n.ens = n.ens)
    if (ncol(time)==0){
      cstv <- function(x,time=time, ...) {
        geoChronR::createSyntheticTimeseries(values = x, ...)
      }
    } else{
      cstv <- function(x,time=time, ...) {
        geoChronR::createSyntheticTimeseries(values = x, time = time[,sample(seq(1,ncol(time)),1)], ...)
      }
    }

    surVals <- purrr::map(valList,
                          cstv,
                          time = time,
                          sameTrend = TRUE,
                          n.ens = ncp)

  }else if(grepl(surrogate.method,pattern = "spectra",ignore.case = T)){

    cstv <- function(x,...){
      g <- is.finite(x)
      out <- surrogateDataFun(ts = x[g],...)
      nc <- ncol(out)
      om <- matrix(NA,nrow = NROW(g),ncol = nc)
      om[g,] <- out
      return(om)
      }

    surVals <- purrr::map(valList,
                          cstv,
                          method = 'ebisuzaki',
                          num_surr = ncp)
  }else if(grepl(surrogate.method,pattern = "shuffle",ignore.case = T)){
    cstv <- function(x,...){
      g <- is.finite(x)
      out <- surrogateDataFun(ts = x[g],...)
      nc <- ncol(out)
      om <- matrix(NA,nrow = NROW(g),ncol = nc)
      om[g,] <- out
      return(om)
    }

    surVals <- purrr::map(valList,
                          cstv,
                          method = 'random_shuffle',
                          num_surr = ncp)
  }

#this way repeats the uncertainty propagation for EACH surrrogate
  pucv <- function(x,...){propagateUncertainty(vals = x,...)}

  if(progress){
    progress <-  glue::glue("Testing null hypothesis with {mc.ens} simulations, each with {n.ens} ensemble members.")
  }

  out <- purrr::map(surVals,
                    pucv,
                    time = time,
                    changeFun = changeFun,
                    n.ens = n.ens,
                    .progress = progress,
                    ...)


#this way uses DIFFERENT surrogates each time through (much faster) But I also don't know that it makes sense
#svm <- purrr::map(surVals,as.data.frame) %>% purrr::list_cbind() %>% as.matrix() %>% suppressMessages()

#out <- propagateUncertainty(time = time, vals = svm,changeFun = changeFun,n.ens = n.ens,...)
#END THIS IDEA

  return(out)

}


simulateAutoCorrelatedUncertainty <- function(sd, n, mean = 0, ar = sqrt(.5),arima.order = c(1,0,0)){
  unc <- arima.sim(list(order = arima.order, ar = ar), n = n)
  unc <- (scale(unc,center = TRUE, scale = TRUE) * sd) + mean
  return(unc)
}
