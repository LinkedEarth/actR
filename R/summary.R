#' Summarize Excursion Output
#'
#' @param object excursion class output
#' @param ... additional inputs
#' @export
summary.excursion <- function(object,...){
  summaryExcursion(object,...)
}

#' Summarize Excursion Output
#'
#' @param object excursion class output
#' @param params.to.print a vector of parameters to print
#' @export
summaryExcursion <- function(object, params.to.print = c("sig.num","n.consecutive","exc.type","min.vals","na.rm")){

  if(!is.na(object$dataSetName) & !is.na(object$paleoData_variableName)){
    title <- glue::glue("{object$dataSetName} - {object$paleoData_variableName}: Excursion test")
  }else if(is.na(object$dataSetName) & !is.na(object$paleoData_variableName)){
    title <- glue::glue("{object$paleoData_variableName}: Excursion test")
  }else if(!is.na(object$dataSetName) & is.na(object$paleoData_variableName)){
    title <- glue::glue("{object$dataSetName}: Excursion test")
  }else{
    title <- glue::glue("Excursion test")
  }



if(object$empirical_pvalue < 0.05){
  resFun <- crayon::green
}else{
  resFun <- crayon::red
}

  ensOut <- object$eventDetection[[1]]
  hasTimeEnsemble <- !identicalVectorsList(ensOut$time)
  hasPaleoEnsemble <- !identicalVectorsList(ensOut$vals)

  if(hasTimeEnsemble){
    hasTimeEnsemble <- crayon::green(hasTimeEnsemble)
  }else{
    hasTimeEnsemble <- crayon::red(hasTimeEnsemble)
  }

  if(hasPaleoEnsemble){
    hasPaleoEnsemble <- crayon::green(hasPaleoEnsemble)
  }else{
    hasPaleoEnsemble <- crayon::red(hasPaleoEnsemble)
  }

if(object$exc.type == "either"){
exc.prefix <- "positive OR negative"
}else if(object$exc.type == "both"){
  exc.prefix <- "positive AND negative"
}else{
  exc.prefix <- object$exc.type
}
  #to do: make "years" terms a variable depending on timeUnits

#print the summary
cat(crayon::bold(glue::glue("{title} results\n\n")))
cat(crayon::silver(glue::glue("Searched for {crayon::bold(exc.prefix)} excursions in a {crayon::bold(object$event.window)} year window around {crayon::bold(object$event.yr)} {object$timeUnits}, with reference windows of {crayon::bold(object$ref.window)} years on either side.\n\n")))
cat("\n")
cat(crayon::bold(glue::glue("Overall result: Empirical p-value = {resFun(object$empirical_pvalue)}\n\n")))
cat(glue::glue("Time uncertainty considered? {hasTimeEnsemble}\n\n"))
cat(glue::glue("Paleo uncertainty considered? {hasPaleoEnsemble}\n\n"))
cat(glue::glue("Error propagation ensemble members = {object$unc.prop.n}\n\n"))
cat(glue::glue("Null hypothesis testing ensemble members = {object$null.hypothesis.n}\n\n"))
cat("\n")
cat(crayon::bold("Parameter choices:\n"))
for(p in params.to.print){
  cat(
    crayon::silver(
      glue::glue(
        "{p} = {object[p]}\n\n"
        )
      )
    )
}



}

#' Summarize shift output
#'
#' @param object shift output
#' @param ... additional inputs (see summaryShift)
#'
#' @return
#' @export
summary.shift <- function(object,...){
  summaryShift(object,...)
}



#' Summarize shift output
#'
#' @param object shift output
#' @param alpha significance level
#' @param params.to.print vector of paramters to print
#'
#' @return
#' @export
summaryShift <- function(object, alpha = 0.05, params.to.print = c("cpt.fun","minimum.segment.length","method","penalty","ncpts.max")){

  params <- createTibbleFromParameterString(object$parameters[1])

  #get shift type
  if(grepl(pattern = "cpt.mean",params$cpt.fun,ignore.case = T) & !grepl(pattern = "cpt.meanVar",params$cpt.fun,ignore.case = T)){
    shift.type <- "Shift in Mean"
  }else if(grepl(pattern = "cpt.var",params$cpt.fun,ignore.case = T)){
    shift.type <- "Shift in Variance"
  }else if(grepl(pattern = "cpt.meanVar",params$cpt.fun,ignore.case = T)){
    shift.type <- "Shift in Mean and Variance"
  }else{
    shift.type <- "Shift"
  }

  if(!is.na(object$input$dataSetName) & !is.na(object$input$paleoData_variableName)){
    title <- glue::glue("{object$input$dataSetName} - {object$input$paleoData_variableName}: {shift.type}")
  }else if(is.na(object$input$dataSetName) & !is.na(object$input$paleoData_variableName)){
    title <- glue::glue("{object$input$paleoData_variableName}: {shift.type}")
  }else if(!is.na(object$input$dataSetName) & is.na(object$input$paleoData_variableName)){
    title <- glue::glue("{object$input$dataSetName}: {shift.type}")
  }else{
    title <- glue::glue("{shift.type}")
  }

  sig.event <- summarizeShiftSignificance(object$shiftDetection,alpha = alpha,paramTib = params)

  n.sig <- nrow(sig.event)

  if(n.sig > 0){
    resFun <- crayon::green
  }else{
    resFun <- crayon::red
  }

  sigMessage <- glue::glue("Detected {resFun(n.sig)} {crayon::italic(shift.type)} event(s) that were significant at the alpha < {alpha} level")


  hasTimeEnsemble <- !identicalMatrixColumns(object$timeEns)
  hasPaleoEnsemble <- !identicalMatrixColumns(object$valEns)

  if(object$time.ens.supplied.n > 1){
    tsmsg <- glue::glue("Time ensemble supplied (n = {object$time.ens.supplied.n})")
  }else{
    tsmsg <- glue::glue("Time ensemble generated in `propagateUncertainties()`)")
  }


  if(object$vals.ens.supplied.n > 1){
    vsmsg <- glue::glue("Values ensemble supplied (n = {object$vals.ens.supplied.n})")
  }else{
    vsmsg <- glue::glue("Values ensemble generated in `propagateUncertainties()`)")
  }

  if(hasTimeEnsemble){
    hasTimeEnsemble <- glue::glue("{crayon::green(hasTimeEnsemble)} {tsmsg}")
  }else{
    hasTimeEnsemble <- crayon::red(hasTimeEnsemble)
  }

  if(hasPaleoEnsemble){
    hasPaleoEnsemble <- glue::glue("{crayon::green(hasPaleoEnsemble)} {vsmsg}")
  }else{
    hasPaleoEnsemble <- crayon::red(hasPaleoEnsemble)
  }

  #to do: make "years" terms a variable depending on timeUnits

  #print the summary
  cat(crayon::bold(glue::glue("{title} results\n\n")))
  cat(crayon::silver(glue::glue("Searched for {shift.type} with a minimum segment length of {crayon::bold(params$minimum.segment.length)} years, summarizing the results over windows of {object$summary.bin.step} years.\n\n")))
  cat("\n")
  cat(crayon::bold(glue::glue("Overall result: {sigMessage}\n\n")))
  if(n.sig > 0 ){
    print(sig.event %>% dplyr::select(time_start,time_end,empirical_pvalue))
  }
  cat(glue::glue("{crayon::bold('Time uncertainty considered?')} {hasTimeEnsemble}\n\n"))
  cat(glue::glue("{crayon::bold('Paleo uncertainty considered?')} {hasPaleoEnsemble}\n\n"))
  cat(glue::glue("Error propagation ensemble members = {object$unc.prop.n}\n\n"))
  cat(glue::glue("Null hypothesis testing ensemble members = {object$null.hypothesis.n} \n  Members simulated using {crayon::italic(object$surrogate.method)} method\n\n"))
  cat("\n")
  cat(crayon::bold("Parameter choices:\n"))
  for(p in params.to.print){
    cat(
      crayon::silver(
        glue::glue(
          "{p} = {params[p]}\n\n"
        )
      )
    )
  }



}
