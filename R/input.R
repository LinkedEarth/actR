#' Prepare input for actR detection methods
#'
#' @param ltt A LiPD-timeseries-tibble, a tibble or data.frame that has the variable(s) of interest, a time variable (age, year or time) along with their metadata, aranged in rows. If ltt = NA, then one in is created from other inputs
#' @param time if ltt is not provided, input a vector or matrix of time (year or age) data. If it's a multicolumn matrix, the columns are time-ensemble members
#' @param vals if ltt is not provided, input a vector or matrix of paleoData. If it's a multicolumn matrix, the columns are value-ensemble members
#' @param time.variable.name If ltt is not provided, specify the name of the time variable (typically 'age' or 'year')
#' @param vals.variable.name If ltt is not provided, specify the name of the paleo variable (e.g., 'd18O' or 'temperature'). Alternatively, if ltt is provided with more rows than expected, this term is used to attempt to select the correct row.
#' @param time.units If ltt is not provided, specify the units the time variable (typically 'yr BP' or 'CE')
#' @param vals.units If ltt is not provided, specify the units the paleo variable (e.g. 'permil' or 'degrees C')
#' @param dataset.name If ltt is not provided, specify the dataset name
#' @param expecting.one.row Are you expecting a one-row tibble? That is, one variable of interest? (default = TRUE)
#' @param sort.by.time Do you want to arrange the output so that the time variable is increasing, while preserving the time-value relationship? (default = TRUE)
#' @param remove.time.nas Do you want to remove observations that are NA, or otherwise not finite, in the time variable (default = TRUE)
#'
#' @return a lipd-ts-tibble ready for analysis in actR
#' @export
prepareInput <- function(ltt = NA,
                         time = NA,
                         vals = NA,
                         time.variable.name = NA,
                         vals.variable.name = NA,
                         time.units = NA,
                         vals.units = NA,
                         dataset.name = NA,
                         expecting.one.row = TRUE,
                         sort.by.time = TRUE,
                         remove.time.nas = TRUE){
  if(!is.na(ltt)){#great, there's already a tibble ts
    #check to make sure it is.
    if(!is.data.frame(ltt)){
      stop(glue::glue("ltt needs to be a data.frame or tibble, but appears to be a {class(ltt)}"))
    }

    nrt <- nrow(ltt)
    if(nrt == 0){stop("This data.frame is empty.")}
    if(expecting.one.row & nrt > 1){#we need to select which row
      if(is.na(vals.variable.name)){
        stop("You entered a multirow tibble, and no vals.variable.name. If you provide vals.variable.name , we'll try to pick the right row; otherwise, supply a one-row tibble of the variable of interest.")
      }
      #try to pick by variableName
      or <- dplyr::filter(ltt,paleoData_variableName == !!vals.variable.name)
      if(nrow(or) == 1){
        print(glue::glue("Selected {vals.variable.name}"))
      }else if(nrow(or) > 1){
        stop(glue::glue("Tried to select row by variable name, but multiple rows matched {vals.variable.name}. You should select the variable of interest and supply a one-row tibble."))
      }else{
        stop(glue::glue("No rows matching {vals.variable.name}"))
      }
      ltt <- or
    }

    #check the contents of the lipd tibble row

    #create time variables if needed
    if(!all(c("time", "timeUnits","timeVariableName") %in% names(ltt))){

      #check age/yearEnsembles
      if("ageensemble" %in% tolower(names(ltt))){
        hasAgeEnsemble <- TRUE
      }else{
        hasAgeEnsemble <- FALSE
      }

      if("yearensemble" %in% tolower(names(ltt))){
        hasYearEnsemble <- TRUE
      }else{
        hasYearEnsemble <- FALSE
      }


      #check age/year
      if("age" %in% names(ltt)){
        hasAge <- TRUE
      }else{
        hasAge <- FALSE
      }

      if("year" %in% names(ltt)){
        hasYear <- TRUE
      }else{
        hasYear <- FALSE
      }

      #prefer age unless told not to
      if(isTRUE(time.variable.name == "year")){
        hasAge <- FALSE
      }
      if(isTRUE(tolower(time.variable.name) == "yearensemble")){
        hasAgeEnsemble <- FALSE
      }

      if(hasAgeEnsemble){
        ltt$time <- ltt$ageEnsemble
        ltt$timeUnits <- ltt$ageUnits
        ltt$timeVariableName <- "age"
      }else if(hasYearEnsemble){
        ltt$time <- ltt$yearEnsemble
        ltt$timeUnits <- ltt$yearUnits
        ltt$timeVariableName <- "year"
      }else if(hasAge){
        ltt$time <- ltt$age
        ltt$timeUnits <- ltt$ageUnits
        ltt$timeVariableName <- "age"
      }else if(hasYear){
        ltt$time <- ltt$year
        ltt$timeUnits <- ltt$yearUnits
        ltt$timeVariableName <- "year"
      }else{
        stop("The LiPD tibble row your provided is missing 'time', 'ageEnsemble', 'yearEnsemble', 'age', and 'year' variables. One of these must be present.")
      }
    }

    if(!"paleoData_values" %in% names(ltt)){
      stop("The LiPD tibble row your provided is missing a 'paleoData_values' variable. This must be included (as a list-column)")
    }
  }else{#looks like we'll be building a barebones lipd tibble ts
    ltt <- tibble::tibble(
      time = list(time),
      timeUnits = time.units,
      timeVariableName = time.variable.name,
      paleoData_values = list(vals),
      paleoData_variableName = vals.variable.name,
      paleoData_units = vals.units,
      dataSetName = dataset.name
    )
    class(ltt) <- c("lipd-ts-tibble",class(ltt))
  }

  #OK, now we have a lipd-ts-tibble, let's check a few things

  #time and paleoData_values
  if(!all(purrr::map_lgl(ltt$time,is.numeric))){
    stop("'time' data need to be a numeric class (matrix is ideal)")
  }
  if(!all(purrr::map_lgl(ltt$paleoData_values,is.numeric))){
    stop("'paleoData_values' data need to be a numeric class (matrix is ideal)")
  }

  #check length
  mf <- function(x,...){
    return(NROW(x$time) == NROW(x$paleoData_values))
  }

  if(!all(apply(ltt,1,mf))){
    stop("'time' and 'paleoData' must have the same number of rows (observations)")
  }

  v2c <- c("timeUnits","timeVariableName","paleoData_units","paleoData_variableName")

  for(v in v2c){
    if(any(is.na(ltt[v]))){
      cat(crayon::magenta(glue::glue("{v} is at least partially absent in the input data\n\n")))
    }
  }


  #sort by time
  for(r in 1:nrow(ltt)){
    if(NCOL(ltt$time[r][[1]]) == 1){#not an ensemble
      if(remove.time.nas){
        wtna <- which(is.finite(ltt$time[r][[1]]))
        ltt$time[r][[1]] <- ltt$time[r][[1]][wtna]
        if(NCOL(ltt$paleoData_values[r][[1]]) == 1){#not an ensemble
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][wtna]
        }else{
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][wtna,]
        }
      }

      if(sort.by.time){
        ss <- sort(ltt$time[r][[1]],index.return = TRUE,decreasing = FALSE)
        ltt$time[r][[1]] <- ltt$time[r][[1]][ss$ix]
        if(NCOL(ltt$paleoData_values[r][[1]]) == 1){#not an ensemble
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][ss$ix]
        }else{
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][ss$ix,]
        }
      }
    }else{#is an ensemble
      medTime <- apply(ltt$time[r][[1]],1,median,na.rm = TRUE)
      if(remove.time.nas){
        wtna <- which(is.finite(medTime))
        ltt$time[r][[1]] <- ltt$time[r][[1]][wtna,]
        if(NCOL(ltt$paleoData_values[r][[1]]) == 1){#not an ensemble
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][wtna]
        }else{
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][wtna,]
        }
      }

      if(sort.by.time){
        ss <- sort(medTime,index.return = TRUE)
        ltt$time[r][[1]] <- ltt$time[r][[1]][ss$ix,]
        if(NCOL(ltt$paleoData_values[r][[1]]) == 1){#not an ensemble
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][ss$ix]
        }else{
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][ss$ix,]
        }
      }

    }
  }

  #check for missing metadata
  if(is.na(ltt$timeUnits)){
    if(NCOL(ltt$time) == 1){
      ltt$timeUnits <- geoChronR::heuristicUnits(ltt$time[[1]])
    }else{
      ltt$timeUnits <- geoChronR::heuristicUnits(ltt$time[[1]][,1])
    }
  }

  if(is.na(ltt$timeVariableName)){
    ltt$timeVariableName <- "time"
  }

  if(is.na(ltt$paleoData_units)){
    ltt$paleoData_units <- "?"
  }
  if(is.na(ltt$paleoData_variableName)){
    ltt$paleoData_variableName <- "unknown"
  }




  return(ltt)
}

