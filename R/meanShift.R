#' Detect a shift in a single time-value pair of paleodata
#'
#' @importFrom digest digest
#' @importFrom changepoint cpt.mean cpt.var cpt.meanvar cpts
#' @importFrom geoChronR gaussianize
#' @param time a vector of time data
#' @param vals a vector paleodata
#' @param minimum.segment.length the minimum allowed length of a detected segment (in time units)
#' @param cpt.fun which function from the changepoint package to use, changepoint::cpt.mean, changepoint::cpt.var or changepoint::cpt.meanvar
#' @param gaussianize Force vals to gaussian distribution before analysis. Default (TRUE). Most (all?) methods in the changepoint package assume gaussian distributions, so this is strongly recommended.
#' @param ... options to pass to cpt.fun . See changepoint function documentation for details.
#'
#' @return A tibble of output data and metadata
#' @export
detectShiftCore = function(time,
                           vals,
                           minimum.segment.length = 1,
                           cpt.fun = changepoint::cpt.mean,
                           gaussianize = TRUE,
                           ...){

  # interpolation options

  #1. interpolate the data to the min resolution of the record
  #TO DO, add more options (med res, binning)
  ad <- abs(diff(time))
  res = min(ad[ad > 0],na.rm = TRUE)
  f = approxfun(time,vals)
  X = seq(min(time,na.rm = TRUE), max(time,na.rm = TRUE), by = res)
  Y = f(X)

  if(gaussianize){
    Y <- as.numeric(geoChronR::gaussianize(Y))
  }

  # pull out ... parameters
  opts <- list(...)

  #check if method can handlle minseglen
  if (isTRUE(opts$method == "SegNeigh") & minimum.segment.length > 1){
    stop("minseglen not yet implemented for SegNeigh method, use PELT instead.")
  }


  #set minimum segment length
  minSeg = ceiling(minimum.segment.length / res)

  # run the mean shift code
  results = cpt.fun(Y,...,minseglen = minSeg)
  resInds = changepoint::cpts(results)


  #get parameters

  params <- glue::glue("cpt.fun = '{deparse(substitute(cpt.fun))}', minimum.segment.length = {minimum.segment.length}, method = '{results@method}', penalty = '{results@pen.type}', pen.value = {results@pen.value}, ncpts.max = {results@ncpts.max}")



  changes <- matrix(X[c(resInds,resInds+1)],ncol = 2,byrow = FALSE)
  if(nrow(changes) > 0){
    change.start <- apply(changes,1,min)
    change.end <- apply(changes,1,max)
  }else{
    change.start <- NA
    change.end <- NA
  }

  # create a hash for each unique time-value pair
  it.hash <- digest::digest(list(time,vals))
  #report out to tibble
  out <- tibble::tibble(time_start = change.start,
                        time_end = change.end,
                        eventDetected = TRUE,
                        eventProbability = NA,
                        time = list(time),
                        vals = list(vals),
                        parameters = as.character(params),
                        it_hash = it.hash) %>%
    new_shiftCore()
  return(out)

} # end function

#
#
# ltt = syntheticTransition
# time = NA
# vals = NA
# time.variable.name = NA
# vals.variable.name = NA
# time.units = NA
# vals.units = NA
# dataset.name = NA
# surrogate.method = "isospectral"
# summary.bin.step = 100
# null.hypothesis.n = 100
# null.quantiles = c(.95,.9)
# time.range = NA

#' Detect a shift in the mean and/or variance of a dataset
#' @description detectShift() allows you to detect a shift in the mean and/or variance of a dataset, and assess its significance given age and data uncertainty relative to a robust null hypothesis. This approach uses the function changepoint::cpt.mean(), changepoint::cpt.var(), or changepoint::cpt.meanvar()  from the changepoint package, propagates inputted or modelled time and/or value ensembles, and summarizes their likelihoods relative to a robust null hypothesis (see ?testNullHypothesis)
#' @inheritParams prepareInput
#' @param summary.bin.step Time interval over which to summarize the results
#' @param null.hypothesis.n How many simulations to run for null hypothesis testing (default = 100)
#' @param null.quantiles What quantiles to report as output from null hypothesis testing (default = c(.95, .9))
#' @inheritParams testNullHypothesis
#' @inheritParams detectShiftCore
#'
#'
#' @return a list that includes error-propagated shift results, null hypothesis testing and metadata
#' @export
detectShift <- function(ltt = NA,
                        time = NA,
                        vals = NA,
                        time.variable.name = NA,
                        vals.variable.name = NA,
                        time.units = NA,
                        vals.units = NA,
                        dataset.name = NA,
                        surrogate.method = "isospectral",
                        summary.bin.step = 100,
                        null.hypothesis.n = 100,
                        null.quantiles = c(.95,.9),
                        time.range = NA,
                        ...){


  #prep inputs
  prepped <- prepareInput(ltt = ltt,
                          time = time,
                          vals = vals,
                          time.variable.name = time.variable.name,
                          vals.variable.name = vals.variable.name,
                          time.units = time.units,
                          vals.units = vals.units,
                          dataset.name = dataset.name,
                          expecting.one.row = TRUE,
                          sort.by.time = TRUE,
                          remove.time.nas = TRUE,
                          time.range = time.range)

  time <- prepped$time[[1]]
  vals <- prepped$paleoData_values[[1]]

  time.ens.supplied.n <- NCOL(time)
  vals.ens.supplied.n <- NCOL(vals)

  #ensemble with uncertainties
  ptm <- proc.time()
  propagated <- propagateUncertainty(time,vals,changeFun = detectShiftCore,...)
  te <- proc.time() - ptm
  te <- te[3]


  propSummary <- summarizeEventProbability(propagated,
                                           bin.step = summary.bin.step,
                                           min.time = min(time, na.rm = T),
                                           max.time = max(time, na.rm = T))
  #null hypothesis
  nh <- testNullHypothesis(time,
                           vals,
                           changeFun = detectShiftCore,
                           surrogate.method = surrogate.method,
                           mc.ens = null.hypothesis.n,
                           how.long.prop = te,
                           ...)


  #get a matrix of nulls
  nhMat <- purrr::map(nh,
                      summarizeEventProbability,
                      bin.step = summary.bin.step,
                      min.time = min(time, na.rm = T),
                      max.time = max(time, na.rm = T)) %>%
    setNames(paste0("nh",seq_along(nh))) %>%
    purrr::map_dfc(purrr::pluck,"event_probability") %>%
    as.matrix()

  nhMat[is.na(nhMat)] <- 0

  #get quantiles for nulls
  nhSummary <- nhMat %>%
    apply(1,quantile,probs = null.quantiles) %>%
    t() %>%
    tibble::as_tibble() %>%
    setNames(paste0("q",null.quantiles))

  #calculate empirical pvalues

  getEmpP <- function(dat,act){
    ef <- ecdf(dat)
    p <- 1 - ef(act)
    p[act == 0] <- 1
    return(p)
  }

  nhSummary$empirical_pvalue <- purrr::array_branch(nhMat,margin = 1) %>%
    purrr::map2_dbl(.y = propSummary$event_probability,.f = getEmpP)

  dsout <- dplyr::bind_cols(propSummary,nhSummary)



  #add in ensemble tables
  ensData <- dplyr::select(propagated,time,vals,it_hash) %>%
    dplyr::group_by(it_hash) %>%
    dplyr::summarize(time = unique(time),
                     vals = unique(vals))

  timeEns <- list2matrix(ensData$time)
  valEns <- list2matrix(ensData$vals)

  #add in metadata
  n.ens<- propagated$nEns[1] #pull example metadata

  out <- list(shiftDetection = dsout,
              parameters = propagated$parameters,
              summary.bin.step = summary.bin.step,
              surrogate.method = surrogate.method,
              unc.prop.n = n.ens,
              null.hypothesis.n = null.hypothesis.n,
              timeEns = timeEns,
              valEns = valEns,
              time.ens.supplied.n = time.ens.supplied.n,
              vals.ens.supplied.n = vals.ens.supplied.n,
              input = prepped) %>%
    new_shift()

  return(out)


}


summarizeShiftSignificance <- function(shiftDetection,alpha = 0.05,paramTib){
  #find significant events
  maxcli <- strsplit(names(shiftDetection),"q") %>%
    purrr::map(pluck,2) %>%
    purrr::map_dbl(~ ifelse(is.null(.x),0,as.numeric(.x))) %>%
    which.max()

  sig.event <- shiftDetection %>%
    dplyr::filter(empirical_pvalue < alpha) %>%
    dplyr::mutate(exceedance = event_probability - .[[maxcli]]) %>%
    dplyr::arrange(empirical_pvalue,dplyr::desc(exceedance))

  #remove those within minimum distance
  bm <- sig.event$time_mid
  if(length(bm) > 1){
    tr <- c()
    for(i in 1:(length(bm) - 1)){
      diffs <- abs(bm-bm[i])
      close <- which(diffs < paramTib$minimum.segment.length)
      ttr <- intersect(close,seq(i+1,length(bm)))
      tr <- c(tr,ttr)
    }

    if(length(tr) > 0){
      sig.event <- sig.event[-tr,]
    }

  }

  return(sig.event)


}
