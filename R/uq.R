
#' Detect a changepoint while propagating proxy and/or age uncertainty
#'
#' @param age
#' @param vals
#' @param changeFun
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
propagateUncertainty <- function(age,
                                 vals,
                                 changeFun,
                                 simulateAgeUncertainty = TRUE,
                                 simulatePaleoUncertainty = TRUE,
                                 n.ens = 100,
                                 bam.model = list(ns = n.ens, name = "bernoulli", param = 0.05),
                                 paleo.uncertainty = sd(vals,na.rm = TRUE),
                                 paleo.ar1 = sqrt(0.5),
                                 paleo.arima.order = c(1,0,0),
                                 summarize = FALSE,

                                 ...){


  #check inputs
  if(!simulateAgeUncertainty & !simulatePaleoUncertainty & NCOL(age) == 1 & NCOL(vals) == 1){
    n.ens <- 1
  }

  nca <- NCOL(age)
  #Prepare age ensemble
  if(nca == 1){#then it's not an ensemble
    #create ensemble?
    if(simulateAgeUncertainty){
      ageMat <- geoChronR::simulateBam(X = matrix(1,nrow = length(age)),
                                       t = as.matrix(age),
                                       model = bam.model,
                                       ageEnsOut = TRUE)$ageEns
    }else{#replicate ages up to n.ens
      ageMat <- matrix(rep(age,n.ens),ncol = n.ens,nrow = NROW(age))
    }
  }else{
    if(nca >= n.ens){
      ageMat <- age[,sample(seq_len(nca),size = n.ens,replace = FALSE)]
    }else if(nca > n.ens){
      ageMat <- age[,sample(seq_len(nca),size = n.ens,replace = TRUE)]
    }
  }

  #make into a list for purrr
  ageList <- purrr::array_branch(ageMat,margin = 2)


  #Now prep paleodata
  ncp <- NCOL(vals)

  if(ncp == 1){#then it's not an ensemble
    #create ensemble?
    if(simulatePaleoUncertainty){
      paleoList <- purrr::rerun(vals + simulateAutoCorrelatedUncertainty(sd = paleo.unc,
                                                                         n = NROW(vals),
                                                                         ar = paleo.ar1,
                                                                         arima.order = paleo.arima.order),.n = n.ens)


    }else{#replicate ages up to n.ens
      paleoList <- matrix(rep(vals,n.ens),ncol = n.ens,nrow = NROW(vals)) %>%
        purrr::array_branch(margin = 2)

    }
  }else{
    if(ncp >= n.ens){
      paleoMat <- vals[,sample(seq_len(ncp),size = n.ens,replace = FALSE)]
    }else if(ncp > n.ens){
      paleoMat <- vals[,sample(seq_len(ncp),size = n.ens,replace = TRUE)]
    }
    paleoList <- purrr::array_branch(paleoMat,margin = 2)
  }


  propagated <- purrr::map2_dfr(ageList,paleoList,changeFun,...)

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
#' @inheritParams detectExcursion
#' @param n.ens The number of ensembles to simulate (default = 100)
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @return a tibble that reports the positivity rate in the synthetics
#' @export
testNullHypothesis <- function(age,
                                  vals,
                                  changeFun,
                                  method = "isospectral",
                                  n.ens = 100,
                                  ...) {

  #create surrogate data for hypothesis testing
  if(grepl(method,pattern = "persis",ignore.case = T)){
    #surVals <- geoChronR::ar1Surrogates(time = age,vals = vals,detrend = TRUE,method = "redfit",n.ens = n.ens)
    surVals <- geoChronR::createSyntheticTimeseries(time = age,values = vals,sameTrend = TRUE,n.ens = n.ens)


  }else if(grepl(method,pattern = "spectra",ignore.case = T)){

  }


  out <- propagateUncertainty(age,surVals,changeFun,n.ens = n.ens,...)

  return(out)

}


simulateAutoCorrelatedUncertainty <- function(sd, n, mean = 0, ar = sqrt(.5),arima.order = c(1,0,0)){
  unc <- arima.sim(list(order = arima.order, ar = ar), n = n)
  unc <- (scale(unc,center = TRUE, scale = TRUE) * sd) + mean
  return(unc)
}
