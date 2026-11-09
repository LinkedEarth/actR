# prepareInput() uses lipdR's time metadata and accepts TS lists

ceRow <- function(ens = FALSE) {
  yrs <- seq(1000, 1990, by = 10)
  r <- tibble::tibble(dataSetName = "ce", paleoData_variableName = "d18O",
                      paleoData_units = "permil",
                      paleoData_values = list(sin(yrs / 40)),
                      year = list(yrs), yearUnits = "CE",
                      time = list(yrs), timeUnits = "yr AD", timeDatum = 0,
                      timeDirection = "prograde", timeExponent = 0)
  if (ens) {
    r$ageEnsemble <- list(matrix(1950 - yrs, nrow = length(yrs), ncol = 4))
    r$ageEnsembleUnits <- "yr BP"
  }
  r
}

test_that("lipdR time is used, and timeVariableName says where it came from", {
  p <- prepareInput(ceRow())
  expect_equal(p$time[[1]], seq(1000, 1990, by = 10))
  expect_equal(p$timeUnits, "yr AD")
  expect_equal(p$timeDirection, "prograde")
  # a hand-built row without lipdR's timeVariableName: time itself is the source
  expect_equal(p$timeVariableName, "time")
})

test_that("time.units converts the time axis", {
  p <- prepareInput(ceRow(), time.units = "yr BP")
  expect_equal(p$time[[1]], sort(1950 - seq(1000, 1990, by = 10)))
  expect_equal(p$timeUnits, "yr BP")
  expect_equal(p$timeDirection, "retrograde")
})

test_that("an age ensemble is preferred, unless a column is named", {
  p <- prepareInput(ceRow(ens = TRUE))
  expect_equal(NCOL(p$time[[1]]), 4)
  expect_equal(p$timeUnits, "yr BP")
  expect_equal(p$timeVariableName, "ageEnsemble")
  p2 <- prepareInput(ceRow(ens = TRUE), time.variable = "year")
  expect_equal(NCOL(p2$time[[1]]), 1)
  expect_equal(p2$timeUnits, "yr AD")
  expect_equal(p2$timeVariableName, "year")
  p3 <- prepareInput(ceRow(ens = TRUE), time.variable = "time")
  expect_equal(NCOL(p3$time[[1]]), 1)
})

test_that("a TS list and an unclassed data frame are accepted", {
  r <- ceRow()
  ts <- purrr::transpose(as.list(r))
  p <- prepareInput(ts)
  expect_equal(p$time[[1]], seq(1000, 1990, by = 10))
  p2 <- prepareInput(as.data.frame(r))
  expect_equal(p2$timeUnits, "yr AD")
})

test_that("old-style input without time metadata still works", {
  r <- tibble::tibble(dataSetName = "old", paleoData_variableName = "x", paleoData_units = "u",
                      paleoData_values = list(1:10 + 0), age = list(seq(10, 100, by = 10)),
                      ageUnits = "yr BP")
  p <- prepareInput(r)
  expect_equal(p$time[[1]], seq(10, 100, by = 10))
  expect_equal(p$timeVariableName, "age")
  r$age <- NULL
  expect_error(prepareInput(r), "time axis")
})

test_that("ka axes are reversed in plots, CE axes are not", {
  expect_true(isRetrogradeAxis("ka"))
  expect_true(isRetrogradeAxis("yr BP"))
  expect_false(isRetrogradeAxis("yr AD"))
  expect_false(isRetrogradeAxis("CE"))
  expect_false(isRetrogradeAxis(NA))
  expect_false(isRetrogradeAxis("?"))
  expect_equal(timeAxisName("yr AD"), "Year")
  expect_equal(timeAxisName("ka"), "Age")
  expect_equal(durationUnits("ka"), "kyr")
  expect_equal(durationUnits("yr BP"), "years")
})

test_that("a missing paleoData_units column is filled in rather than an error", {
  r <- ceRow()
  r$paleoData_units <- NULL
  expect_equal(prepareInput(r)$paleoData_units, "?")
})
