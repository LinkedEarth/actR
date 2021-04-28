detectShiftCore = function(age,
                           vals,
                           minimum.segment.length = 1,
                           cpt.fun = changepoint::cpt.mean,
                           ...){


  ## Written by Hannah Kolus, 09/04/2018
  ## Searches for mean shifts in the specified record, using the changepoint package function cpt.mean().
  ## In addition, tests for significant differences in means from the results of cpt.mean(), with the
  ## potential to reject some changepoints as insignificant.

  # interpolation options

  #1. interpolate the data to the min resolution of the record
  #TO DO, add more options (med res, binning)
  ad <- abs(diff(age))
  res = min(ad[ad > 0],na.rm = TRUE)
  f = approxfun(age,vals)
  X = seq(min(age,na.rm = TRUE), max(age,na.rm = TRUE), by = res)
  Y = f(X)


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
  change.start <- apply(changes,1,min)
  change.end <- apply(changes,1,max)

  it.hash <- digest::digest(list(age,vals))
  #report out to tibble
  out <- tibble::tibble(time_start = change.start,
                        time_end = change.end,
                        eventDetected = TRUE,
                        eventProbability = NA,
                        age = list(age),
                        vals = list(vals),
                        parameters = as.character(params),
                        it_hash = it.hash)
  return(out)

} # end function

detectShift <- function(age,
                        vals,
                        surrogate.method = "isospectral",
                        summary.bin.step = 100,
                        null.hypothesis.n = 100,
                        null.quantiles = c(.9,.95),
                        ...){

  #ensemble with uncertainties
  ptm <- proc.time()
  propagated <- propagateUncertainty(age,vals,changeFun = detectShiftCore,...)
  te <- proc.time() - ptm
  te <- te[3]


  propSummary <- summarizeEventProbability(propagated,
                                           bin.step = summary.bin.step,
                                           min.age = min(age, na.rm = T),
                                           max.age = max(age, na.rm = T))
  #null hypothesis
  nh <- testNullHypothesis(age,
                           vals,
                           changeFun = detectShiftCore,
                           surrogate.method = surrogate.method,
                           mc.ens = null.hypothesis.n,
                           how.long.prop = te,
                           ...)


  nhSummary <- purrr::map(nh,
                          summarizeEventProbability,
                          bin.step = summary.bin.step,
                          min.age = min(age, na.rm = T),
                          max.age = max(age, na.rm = T)) %>%
    setNames(paste0("nh",seq_along(nh))) %>%
    purrr::map_dfc(purrr::pluck,"event_probability") %>%
    as.matrix() %>%
    apply(1,quantile,probs = null.quantiles) %>%
    t() %>%
    tibble::as_tibble() %>%
    setNames(paste0("q",null.quantiles))

  dsout <- dplyr::bind_cols(propSummary,nhSummary)



  #add in ensemble tables
  ensData <- dplyr::select(propagated,age,vals,it_hash) %>%
    dplyr::group_by(it_hash) %>%
    dplyr::summarize(age = unique(age),
                     vals = unique(vals))

  ageEns <- list2matrix(ensData$age)
  valEns <- list2matrix(ensData$vals)

  #add in metadata

  out <- list(shiftDetection = dsout,
              parameters = propagated$parameters,
              null.hypothesis.n = null.hypothesis.n,
              ageEns = ageEns,
              valEns = valEns)

return(out)


}
