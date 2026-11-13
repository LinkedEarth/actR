#' Prepare input for actR detection methods
#'
#' @param ltt A LiPD-timeseries-tibble (`lipdR::ts2tibble()`), a TS list (`lipdR::extractTs()`), or a tibble or data.frame that has the variable(s) of interest, a time variable (age, year or time) along with their metadata, arranged in rows. If ltt = NA, then one in is created from other inputs. The time axis is an age or year ensemble when there is one, then lipdR's `time` (the primary time column), then `age`, then `year`.
#' @param time if ltt is not provided, input a vector or matrix of time (year or age) data. If it's a multicolumn matrix, the columns are time-ensemble members
#' @param vals if ltt is not provided, input a vector or matrix of paleoData. If it's a multicolumn matrix, the columns are value-ensemble members
#' @param time.variable Which time variable to use. If ltt is provided, the name of a time column (e.g. 'age', 'year', 'ageEnsemble', or any column lipdR found to be a time column; see `lipdR::tsTimeVariables()`), passed to `lipdR::getTsTime(variable =)`. NA (default) uses an age or year ensemble when there is one, then lipdR's `time`, the record's primary time axis. If ltt is not provided, a name for the time variable (typically 'age' or 'year').
#' @param time.variable.name Deprecated; use `time.variable`.
#' @param vals.variable.name If ltt is not provided, specify the name of the paleo variable (e.g., 'd18O' or 'temperature'). Alternatively, if ltt is provided with more rows than expected, this term is used to attempt to select the correct row.
#' @param time.units If ltt is not provided, the units of `time` (e.g. 'yr BP' or 'CE'). If ltt is provided, the time scale (e.g. 'yr BP', 'ka' or 'CE') to put the time axis on with `lipdR::getTsTime()`, so that event years, windows and bins are on that scale whatever the record's own units. NA (default) keeps the record's own time scale.
#' @param vals.units If ltt is not provided, specify the units the paleo variable (e.g. 'permil' or 'degrees C')
#' @param dataset.name If ltt is not provided, specify the dataset name
#' @param expecting.one.row Are you expecting a one-row tibble? That is, one variable of interest? (default = TRUE)
#' @param sort.by.time Do you want to arrange the output so that the time variable is increasing, while preserving the time-value relationship? (default = TRUE)
#' @param remove.time.nas Do you want to remove observations that are NA, or otherwise not finite, in the time variable (default = TRUE)
#' @param time.range Optionally enter a time range (as minimum and maximum) that you'd like to restrict the analysis to. (default = NA)
#' @param warn.assumed Warn when a record's time scale was assumed, because its units don't say? (default = TRUE)
#'
#' @return a lipd-ts-tibble ready for analysis in actR
#' @export
prepareInput <- function(ltt = NA,
                         time = NA,
                         vals = NA,
                         time.variable = NA,
                         vals.variable.name = NA,
                         time.units = NA,
                         vals.units = NA,
                         dataset.name = NA,
                         expecting.one.row = TRUE,
                         sort.by.time = TRUE,
                         remove.time.nas = TRUE,
                         time.range = NA,
                         warn.assumed = TRUE,
                         time.variable.name = NA){
  time.variable <- resolveTimeVariable(time.variable, time.variable.name)
  if(!is.data.frame(ltt) && is.list(ltt)){#a TS list (lipdR::extractTs()), so make it a tibble
    ltt <- lipdR::as.lipdTsTibble(ltt)
  }
  if(!all(is.na(ltt))){#great, there's already a tibble ts
    #check to make sure it is.
    if(!is.data.frame(ltt)){#change this to class once implemented in lipdR
      stop(glue::glue("ltt needs to be a data.frame or tibble, but appears to be a {class(ltt)}"))
    }

    nrt <- nrow(ltt)
    if(nrt == 0){stop("This data.frame is empty.")}
    if(expecting.one.row & nrt > 1){#we need to select which row
      if(is.na(vals.variable.name)){
        cat(crayon::blue("Which variable do you want to test? (Specify this in 'vals.variable.name' to avoid this message)\n\n"))
        for(p in 1:nrt){
          cat(glue::glue("{p} - {ltt$paleoData_variableName[p]} ({n = length(ltt$paleoData_values[[p]])})\n\n"))
        }
        n <- as.numeric(readline(prompt="Please type the number for the correct match: "))
        ltt <- ltt[n,]

      }else{
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
    }

    #check the contents of the lipd tibble row

    #create the time axis from lipdR's time metadata (or age/year),
    #converting to time.units if given
    ltt <- addActTime(ltt, time.variable = time.variable, time.units = time.units,
                      warn.assumed = warn.assumed)

    if(!"paleoData_values" %in% names(ltt)){
      stop("The LiPD tibble row your provided is missing a 'paleoData_values' variable. This must be included (as a list-column)")
    }
  }else{#looks like we'll be building a barebones lipd tibble ts
    ltt <- tibble::tibble(
      time = list(time),
      timeUnits = time.units,
      timeVariableName = time.variable,
      timeActR = TRUE,
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
    if(!v %in% names(ltt)){
      ltt[[v]] <- NA_character_
    }
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

  #truncate by time
  if(!any(is.na(time.range))){
    for(r in 1:nrow(ltt)){
      if(NCOL(ltt$time[r][[1]]) == 1){#not an ensemble
        good.time <- which(ltt$time[r][[1]] >= min(time.range) & ltt$time[r][[1]] <= max(time.range))
        if(length(good.time) < 1){stop("No values left after truncating for time.range")}
        ltt$time[r][[1]] <- ltt$time[r][[1]][good.time]
        if(NCOL(ltt$paleoData_values[r][[1]]) == 1){#not an ensemble
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][good.time]
        }else{
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][good.time,]
        }
      }else{#is an ensemble
        medTime <- apply(ltt$time[r][[1]],1,median,na.rm = TRUE)
        good.time <- which(medTime >= min(time.range) & medTime <= max(time.range))
        if(length(good.time) < 1){stop("No values left after truncating for time.range")}
        ltt$time[r][[1]] <- ltt$time[r][[1]][good.time,]
        if(NCOL(ltt$paleoData_values[r][[1]]) == 1){#not an ensemble
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][good.time]
        }else{
          ltt$paleoData_values[r][[1]] <- ltt$paleoData_values[r][[1]][good.time,]
        }
      }
    }

    #truncate the original age data too, to avoid issues
    for(r in 1:nrow(ltt)){

      if(length(ltt$ageEnsemble[[r]]) > 0){
        ltt$ageEnsemble[[r]]<- ltt$ageEnsemble[[r]][good.time,]
      }
      if(length(ltt$yearEnsemble[[r]]) > 0){
        ltt$yearEnsemble[[r]] <- ltt$yearEnsemble[[r]][good.time,]
      }
      if(length(ltt$age[[r]]) > 0){
        ltt$age[[r]] <- ltt$age[[r]][good.time]
      }
      if(length(ltt$year[[r]]) > 0){
        ltt$year[[r]] <- ltt$year[[r]][good.time]
      }
    }
  }

  #check for missing metadata
  for(r in 1:nrow(ltt)){
    if(is.na(ltt$timeUnits[r])){
      tm <- if(NCOL(ltt$time[[r]]) == 1) ltt$time[[r]] else ltt$time[[r]][,1]
      ltt$timeUnits[r] <- ens::heuristicUnits(as.numeric(tm))
      if(warn.assumed){
        if(isTRUE(grepl("guess", ltt$timeUnits[r]))){
          warning("No time units were given, and they couldn't be guessed from the values. Set time.units to say.",
                  call. = FALSE)
        }else{
          warning("No time units were given, so ", ltt$timeUnits[r],
                  " was guessed from the values. Set time.units to say.", call. = FALSE)
        }
      }
    }
  }

  if(any(is.na(ltt$timeVariableName))){
    ltt$timeVariableName[is.na(ltt$timeVariableName)] <- "time"
  }

  if(any(is.na(ltt$paleoData_units))){
    ltt$paleoData_units[is.na(ltt$paleoData_units)] <- "?"
  }
  if(any(is.na(ltt$paleoData_variableName))){
    ltt$paleoData_variableName[is.na(ltt$paleoData_variableName)] <- "unknown"
  }




  return(ltt)
}




#' Add actR time columns to a timeseries tibble
#'
#' Uses `lipdR::getTsTime()` on each row, so the time axis is read from the
#' record's time columns and their time scales, converted to `time.units`. Rows
#' that actR has already prepared (column `timeActR`) keep their `time`, which
#' is only converted, because it has been sorted and trimmed with the values.
#'
#' @param ltt a timeseries tibble
#' @param time.variable optional time variable to use (see prepareInput())
#' @param time.units optional time scale to convert to (NA keeps each record's own scale)
#' @param warn.assumed give one warning for records whose time scale was assumed?
#' @return ltt with time, timeUnits, timeDatum, timeDirection, timeExponent,
#'   timeVariableName (the source of the time axis) and timeActR columns
#' @keywords internal
addActTime <- function(ltt, time.variable = NA, time.units = NA, warn.assumed = TRUE){
  rows <- purrr::transpose(as.list(ltt))
  tv <- if(is.na(time.variable)) NA_character_ else as.character(time.variable)
  prepared <- "timeActR" %in% names(ltt)
  assumed <- character(0)
  one <- function(r){
    label <- paste(stats::na.omit(c(unlist(r$dataSetName)[1], unlist(r$paleoData_variableName)[1])), collapse = " ")
    if(prepared && isTRUE(r$timeActR)){
      #already an actR input: keep its time axis (no pointer), only converting
      r <- r[intersect(names(r), c("time", "timeUnits", "timeDatum", "timeDirection", "timeExponent", "timeSource", "timeVariableName"))]
      src <- r$timeVariableName
      r$timeVariableName <- NULL
      if(is.null(r$timeUnits) || is.na(r$timeUnits)){
        if(!is.na(time.units)){
          stop("The time units of ", if(nzchar(label)) label else "this timeseries",
               " are unknown, so it can't be put on ", time.units, ".", call. = FALSE)
        }
        tm <- r$time
        attr(tm, "timeUnits") <- NA_character_
        attr(tm, "timeVariableName") <- src %||% "time"
        return(tm)
      }
      tm <- lipdR::getTsTime(r, units = time.units, ensemble = FALSE, warn.assumed = FALSE)
      attr(tm, "timeVariableName") <- src %||% "time"
      return(tm)
    }
    tm <- tryCatch({
      if(is.na(tv)){
        lipdR::getTsTime(r, units = time.units, ensemble = TRUE, warn.assumed = FALSE)
      }else{
        lipdR::getTsTime(r, units = time.units, ensemble = FALSE, variable = tv, warn.assumed = FALSE)
      }
    }, error = function(e){
      if(!is.na(tv) && grepl("not a time variable", conditionMessage(e)) && !is.null(r[[tv]])){
        #a column lipdR doesn't know as time: use it with its units, or time.units
        out <- lipdR::getTsTime(list(time = r[[tv]], timeUnits = r[[paste0(tv, "Units")]] %||% time.units),
                                units = time.units, ensemble = FALSE, warn.assumed = FALSE)
        attr(out, "timeVariableName") <- tv
        return(out)
      }
      stop("Can't find the time axis of ", if(nzchar(label)) label else "this timeseries", ": ", conditionMessage(e),
           " ('time', 'ageEnsemble', 'yearEnsemble', 'age' or 'year' must be present.)", call. = FALSE)
    })
    if(identical(attr(tm, "timeSource"), "assumed")){
      assumed <<- c(assumed, if(nzchar(label)) label else "a timeseries")
    }
    tm
  }
  out <- purrr::map(rows, one)
  if(warn.assumed && length(assumed) > 0){
    warning("The time scale of ", length(assumed), " timeseries was assumed, because their units don't say: ",
            paste(unique(assumed), collapse = "; "),
            ". Set their units (or timeDatum, timeDirection and timeExponent) to say.", call. = FALSE)
  }
  units <- purrr::map_chr(out, \(x) attr(x, "timeUnits") %||% NA_character_)
  ltt$timeVariableName <- purrr::map_chr(out, \(x) attr(x, "timeVariableName") %||% "time")
  ltt$time <- purrr::map(out, \(x){
    attributes(x) <- attributes(x)[intersect(names(attributes(x)), "dim")]
    if(NCOL(x) == 1) as.numeric(x) else x
  })
  ltt$timeUnits <- units
  scales <- purrr::map(units, \(u) if(is.na(u)) NULL else lipdR::timeScale(u))
  ltt$timeDatum <- purrr::map_dbl(scales, \(x) x$timeDatum %||% NA_real_)
  ltt$timeDirection <- purrr::map_chr(scales, \(x) x$timeDirection %||% NA_character_)
  ltt$timeExponent <- purrr::map_dbl(scales, \(x) x$timeExponent %||% NA_real_)
  ltt$timeActR <- TRUE
  ltt
}

#' Resolve time.variable and its deprecated alias time.variable.name
#' @param time.variable the time.variable argument
#' @param time.variable.name the deprecated time.variable.name argument
#' @return the time variable to use
#' @keywords internal
resolveTimeVariable <- function(time.variable, time.variable.name){
  if(!all(is.na(time.variable.name))){
    warning("time.variable.name is deprecated; use time.variable.", call. = FALSE)
    if(all(is.na(time.variable))){
      time.variable <- time.variable.name
    }
  }
  time.variable
}

`%||%` <- function(x, y){
  if(is.null(x)) y else x
}

#' Is a time axis retrograde (counting back, like yr BP or ka), so plots run right to left?
#' @param units a time units string
#' @return logical; FALSE for missing or unknown ("?") units
#' @keywords internal
isRetrogradeAxis <- function(units){
  u <- unlist(units)[1]
  if(is.null(u) || is.na(u) || !nzchar(u) || u == "?"){
    return(FALSE)
  }
  !lipdR::isPrograde(u)
}

#' Is a time axis prograde (counting forward, like year CE)?
#' @inheritParams isRetrogradeAxis
#' @return logical; FALSE for missing or unknown units
#' @keywords internal
isProgradeAxis <- function(units){
  u <- unlist(units)[1]
  if(is.null(u) || is.na(u) || !nzchar(u) || u == "?"){
    return(FALSE)
  }
  lipdR::isPrograde(u)
}

#' Name of the time axis, from its direction: "Year" counts forward, "Age" back
#' @inheritParams isRetrogradeAxis
#' @keywords internal
timeAxisName <- function(units){
  if(isProgradeAxis(units)) "Year" else if(isRetrogradeAxis(units)) "Age" else "Time"
}

#' Units of a duration on a time axis (years, kyr, Myr)
#' @inheritParams isRetrogradeAxis
#' @keywords internal
durationUnits <- function(units){
  u <- unlist(units)[1]
  e <- if(is.null(u) || is.na(u) || u == "?") 0 else tryCatch(lipdR::timeScale(u)$timeExponent, error = function(e) 0)
  switch(as.character(e), "3" = "kyr", "6" = "Myr", "9" = "Gyr", "years")
}
