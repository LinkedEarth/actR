# actR works on any time axis: CE and BP records give the same answer, in the
# right direction, and the time scale is chosen or defaulted consistently.

# a record whose values step up by 3 at 1500 CE (450 yr BP): younger is higher
shiftRecord <- function(axis = c("CE", "BP")){
  axis <- match.arg(axis)
  set.seed(7)
  yrs <- seq(1000, 1990, by = 10)
  vals <- rnorm(length(yrs), sd = 0.5) + ifelse(yrs > 1500, 3, 0)
  if(axis == "CE"){
    tibble::tibble(dataSetName = "ce", paleoData_variableName = "x", paleoData_units = "u",
                   paleoData_values = list(vals), year = list(yrs), yearUnits = "CE",
                   time = list(yrs), timeUnits = "yr AD", timeDatum = 0,
                   timeDirection = "prograde", timeExponent = 0)
  }else{
    tibble::tibble(dataSetName = "bp", paleoData_variableName = "x", paleoData_units = "u",
                   paleoData_values = list(vals), age = list(1950 - yrs), ageUnits = "yr BP",
                   time = list(1950 - yrs), timeUnits = "yr BP", timeDatum = 1950,
                   timeDirection = "retrograde", timeExponent = 0)
  }
}

runShift <- function(r, ...){
  set.seed(11)
  suppressMessages(capture.output(out <- detectShift(r, n.ens = 10, null.hypothesis.n = 10,
                                                    summary.bin.step = 100,
                                                    minimum.segment.length = 100, ...)))
  out
}

peak <- function(out) out[which.max(out$event_probability), ]

test_that("a CE and a BP record give the same shift, in the right direction", {
  ce <- runShift(shiftRecord("CE"))
  bp <- runShift(shiftRecord("BP"))

  expect_equal(unique(ce$timeUnits), "yr AD")
  expect_equal(unique(bp$timeUnits), "yr BP")
  expect_equal(unique(ce$time_direction), "prograde")
  expect_equal(unique(bp$time_direction), "retrograde")

  # the shift is found at 1500 CE = 450 yr BP
  expect_true(abs(peak(ce)$time_mid - 1500) <= 100)
  expect_true(abs(peak(bp)$time_mid - 450) <= 100)

  # values went up, forward in time: a positive shift, with a positive delta, on both axes
  expect_gt(peak(ce)$event_probability_positive, peak(ce)$event_probability_negative)
  expect_gt(peak(bp)$event_probability_positive, peak(bp)$event_probability_negative)
  expect_gt(peak(ce)$deltas, 0)
  expect_gt(peak(bp)$deltas, 0)
  expect_equal(peak(ce)$deltas, peak(bp)$deltas, tolerance = 0.1)

  # summary window widths are positive bin widths, on both axes
  expect_equal(unique(ce$summary.bin.step), 100)
  expect_equal(unique(bp$summary.bin.step), 100)
  expect_output(summary(ce), "windows of 100 years \\(yr AD\\)")
})

test_that("time.units puts a CE record on a BP axis, with the same answer", {
  ce2bp <- runShift(shiftRecord("CE"), time.units = "yr BP")
  expect_equal(unique(ce2bp$timeUnits), "yr BP")
  expect_true(abs(peak(ce2bp)$time_mid - 450) <= 100)
  expect_gt(peak(ce2bp)$event_probability_positive, peak(ce2bp)$event_probability_negative)
  expect_gt(peak(ce2bp)$deltas, 0)
})

test_that("time.variable picks a time column other than the primary one", {
  r <- shiftRecord("CE")
  # as lipdR >= 0.8.0 extracts it: time points to year, and age carries its own scale
  r$timeVariableName <- "year"
  r$yearTimeUnits <- "yr AD"; r$yearTimeDatum <- 0; r$yearTimeDirection <- "prograde"
  r$yearTimeExponent <- 0; r$yearTimeSource <- "units"
  r$age <- list(rev(seq(0.04, 1.03, by = 0.01)))
  r$ageUnits <- "ka"
  r$ageTimeUnits <- "ka"; r$ageTimeDatum <- 1950; r$ageTimeDirection <- "retrograde"
  r$ageTimeExponent <- 3; r$ageTimeSource <- "units"

  p <- prepareInput(r)
  expect_equal(p$timeVariableName, "year")
  expect_equal(p$timeUnits, "yr AD")

  pa <- prepareInput(r, time.variable = "age")
  expect_equal(pa$timeVariableName, "age")
  expect_equal(pa$timeUnits, "ka")
  expect_equal(pa$time[[1]], seq(0.04, 1.03, by = 0.01))

  expect_warning(pd <- prepareInput(r, time.variable.name = "age"), "deprecated")
  expect_equal(pd$timeVariableName, "age")
})

test_that("records on different scales are put on yr BP, with a message; a shared scale is kept", {
  two <- dplyr::bind_rows(shiftRecord("CE"), shiftRecord("BP"))
  expect_message(
    out <- suppressWarnings(detectMultipleExcursions(two, n.ens = 5, null.hypothesis.n = 5,
                                                     event.yr = 450, event.window = 100,
                                                     ref.window = 200, seed = 1)),
    "different time scales")
  expect_true(all(out$timeUnits == "yr BP"))

  ce2 <- dplyr::bind_rows(shiftRecord("CE"), shiftRecord("CE"))
  ce2$dataSetName <- c("ce1", "ce2")
  expect_no_message(
    out2 <- suppressWarnings(detectMultipleExcursions(ce2, n.ens = 5, null.hypothesis.n = 5,
                                                      event.yr = 1500, event.window = 100,
                                                      ref.window = 200, seed = 1)),
    message = "different time scales")
  expect_true(all(out2$timeUnits == "yr AD"))
})

test_that("assumed time scales give one warning per call", {
  noUnits <- function(name){
    r <- shiftRecord("BP")
    r$dataSetName <- name
    r$ageUnits <- NULL
    r[c("time", "timeUnits", "timeDatum", "timeDirection", "timeExponent")] <- NULL
    r
  }
  two <- dplyr::bind_rows(noUnits("a"), noUnits("b"))
  w <- character(0)
  withCallingHandlers(prepareInput(two, expecting.one.row = FALSE),
                      warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  assumedW <- grep("assumed", w, value = TRUE)
  expect_length(assumedW, 1)
  expect_match(assumedW, "2 timeseries")
  expect_no_warning(prepareInput(two, expecting.one.row = FALSE, warn.assumed = FALSE))

  w2 <- character(0)
  withCallingHandlers(
    suppressMessages(detectMultipleExcursions(two, n.ens = 5, null.hypothesis.n = 5,
                                              event.yr = 450, event.window = 100,
                                              ref.window = 200, seed = 1)),
    warning = function(cnd){ w2 <<- c(w2, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_length(grep("assumed", w2), 1)
})
