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

test_that("lipdR time is used, with a timeVariableName from its direction", {
  p <- prepareInput(ceRow())
  expect_equal(p$time[[1]], seq(1000, 1990, by = 10))
  expect_equal(p$timeUnits, "yr AD")
  expect_equal(p$timeVariableName, "year")
})

test_that("time.units converts the time axis", {
  p <- prepareInput(ceRow(), time.units = "yr BP")
  expect_equal(p$time[[1]], sort(1950 - seq(1000, 1990, by = 10)))
  expect_equal(p$timeUnits, "yr BP")
  expect_equal(p$timeVariableName, "age")
})

test_that("an age ensemble is preferred, unless a column is named", {
  p <- prepareInput(ceRow(ens = TRUE))
  expect_equal(NCOL(p$time[[1]]), 4)
  expect_equal(p$timeUnits, "yr BP")
  p2 <- prepareInput(ceRow(ens = TRUE), time.variable.name = "year")
  expect_equal(NCOL(p2$time[[1]]), 1)
  expect_equal(p2$timeUnits, "yr AD")
  p3 <- prepareInput(ceRow(ens = TRUE), time.variable.name = "time")
  expect_equal(NCOL(p3$time[[1]]), 1)
})

test_that("a TS list and an unclassed data frame are accepted", {
  r <- ceRow()
  ts <- purrr::transpose(as.list(r))
  p <- prepareInput(ts)
  expect_equal(p$time[[1]], seq(1000, 1990, by = 10))
  p2 <- prepareInput(as.data.frame(r))
  expect_equal(p2$timeVariableName, "year")
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

test_that("ka axes are reversed in plots", {
  expect_true(isRetrogradeUnits("ka"))
  expect_true(isRetrogradeUnits("yr BP"))
  expect_false(isRetrogradeUnits("yr AD"))
  expect_false(isRetrogradeUnits("CE"))
})

test_that("a missing paleoData_units column is filled in rather than an error", {
  r <- ceRow()
  r$paleoData_units <- NULL
  expect_equal(prepareInput(r)$paleoData_units, "?")
})
