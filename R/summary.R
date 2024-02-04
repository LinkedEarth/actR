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


  if(!is.na(object$pvalue)){
    if(object$pvalue < 0.05){
      resFun <- crayon::green
    }else{
      resFun <- crayon::red
    }
  }else{
    resFun <- crayon::red
    object$pvalue <- "Excursion test failed to run."
  }

  ensOut <- object$event_detection[[1]]
  hasTimeEnsemble <- !identicalVectorsList(ensOut$time)
  hasPaleoEnsemble <- !identicalVectorsList(ensOut$vals)
  hasParamEnsemble <-  hasParameterEnsemble(object)

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

  if(hasParamEnsemble){
    hasParamEnsemblePrint <- crayon::green(hasParamEnsemble)



  }else{
    hasParamEnsemblePrint <- crayon::red(hasParamEnsemble)
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
cat(crayon::bold(glue::glue("Overall result: Empirical p-value = {resFun(round(object$pvalue,3))}\n\n")))
cat(glue::glue("Time uncertainty considered? {hasTimeEnsemble}\n\n"))
cat(glue::glue("Paleo uncertainty considered? {hasPaleoEnsemble}\n\n"))
if(hasParamEnsemble){
  cat(glue::glue("Parametric uncertainty considered? {hasParamEnsemblePrint}: {crayon::silver(crayon::italic(whichParametersInEnsemble(object)))}\n\n"))
}else{
  cat(glue::glue("Parametric uncertainty considered? {hasParamEnsemblePrint}\n\n"))
}
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
#' @export
summary.shift <- function(object,...){
  summaryShift(object,...)
}



#' Summarize shift output
#'
#' @param object shift output
#' @param alpha significance level
#' @param params.to.print vector of paramters to print
#' @export
summaryShift <- function(object,
                         alpha = 0.05,
                         params.to.print = c("cpt.fun","minimum.segment.length","method","penalty","ncpts.max")){

  #params <- createTibbleFromParameterString(object$parameters[1])

  #get shift type
  if(grepl(pattern = "cpt.mean",object$cpt.fun[1],ignore.case = T) & !grepl(pattern = "cpt.meanVar",object$cpt.fun[1],ignore.case = T)){
    shift.type <- "Shift in Mean"
  }else if(grepl(pattern = "cpt.var",object$cpt.fun[1],ignore.case = T)){
    shift.type <- "Shift in Variance"
  }else if(grepl(pattern = "cpt.meanVar",object$cpt.fun[1],ignore.case = T)){
    shift.type <- "Shift in Mean and Variance"
  }else{
    shift.type <- "Shift"
  }

  if(!is.na(object$dataSetName[1]) & !is.na(object$paleoData_variableName[1])){
    title <- glue::glue("{object$dataSetName[1]} - {object$paleoData_variableName[1]}: {shift.type}")
  }else if(is.na(object$input$dataSetName) & !is.na(object$input$paleoData_variableName)){
    title <- glue::glue("{object$paleoData_variableName[1]}: {shift.type}")
  }else if(!is.na(object$input$dataSetName) & is.na(object$input$paleoData_variableName)){
    title <- glue::glue("{object$dataSetName[1]}: {shift.type}")
  }else{
    title <- glue::glue("{shift.type}")
  }

  sig.event <- summarizeShiftSignificance(object,alpha = alpha,minimum.segment.length = object$minimum.segment.length[1])

  n.sig <- nrow(sig.event)

  if(n.sig > 0){
    resFun <- crayon::green
  }else{
    resFun <- crayon::red
  }

  sigMessage <- glue::glue("Detected {resFun(n.sig)} {crayon::italic(shift.type)} event(s) that were significant at the alpha < {alpha} level")


  hasTimeEnsemble <- !identicalMatrixColumns(object$time[[1]])
  hasPaleoEnsemble <- !identicalMatrixColumns(object$paleoData_values[[1]])

  if(object$time.ens.supplied.n[1] > 1){
    tsmsg <- glue::glue("Time ensemble supplied (n = {object$time.ens.supplied.n[1]})")
  }else{
    tsmsg <- glue::glue("Time ensemble generated in `propagateUncertainties()`)")
  }


  if(object$vals.ens.supplied.n[1] > 1){
    vsmsg <- glue::glue("Values ensemble supplied (n = {object$vals.ens.supplied.n[1]})")
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
  cat(crayon::silver(glue::glue("Searched for {shift.type} with a minimum segment length of {crayon::bold(object$minimum.segment.length[1])} years, summarizing the results over windows of {object$summary.bin.step[1]} years.\n\n")))
  cat("\n")
  cat(crayon::bold(glue::glue("Overall result: {sigMessage}\n\n")))
  if(n.sig > 0 ){
    print(sig.event %>% dplyr::select(time_start,time_end,pvalue,deltas),n=nrow(sig.event))
  }
  cat(glue::glue("{crayon::bold('Time uncertainty considered?')} {hasTimeEnsemble}\n\n"))
  cat(glue::glue("{crayon::bold('Paleo uncertainty considered?')} {hasPaleoEnsemble}\n\n"))
  cat(glue::glue("Error propagation ensemble members = {object$unc.prop.n[1]}\n\n"))
  cat(glue::glue("Null hypothesis testing ensemble members = {object$null.hypothesis.n[1]} \n  Members simulated using {crayon::italic(object$surrogate.method[1])} method\n\n"))
  cat("\n")
  cat(crayon::bold("Parameter choices:\n"))
  for(p in params.to.print){
    cat(
      crayon::silver(
        glue::glue(
          "{p} = {object[1,p]}\n\n"
        )
      )
    )
  }



}
