# Golden-reference characterization tests. Reference outputs captured by
# tools/capture_reference.R from actR v0.2.2 (pre-stack-rewiring). These lock
# the numerical behavior of the core detectors, uncertainty propagation, and
# null-hypothesis testing through the migration to the ens/lipdViz stack.

ref <- readRDS(test_path("reference_outputs.rds"))

makeSyntheticExcursion <- function() {
  set.seed(42)
  time <- seq(0, 1000, by = 10)
  vals <- rnorm(length(time))
  vals[48:52] <- vals[48:52] + 4
  list(time = time, vals = vals)
}

makeSyntheticShift <- function() {
  set.seed(43)
  time <- seq(0, 1000, by = 10)
  vals <- rnorm(length(time)) + ifelse(time > 500, 2, 0)
  list(time = time, vals = vals)
}

test_that("detectExcursionCore matches pre-migration reference", {
  exc <- makeSyntheticExcursion()
  out <- detectExcursionCore(exc$time, exc$vals,
                             event.yr = 500, event.window = 100,
                             ref.window = 200)
  expect_equal(out, ref$excursionCore)
  expect_true(out$eventDetected)
})

test_that("detectShiftCore matches pre-migration reference", {
  shf <- makeSyntheticShift()
  out <- detectShiftCore(shf$time, shf$vals, minimum.segment.length = 100)
  expect_equal(out, ref$shiftCore)
})

test_that("propagateUncertainty matches pre-migration reference", {
  exc <- makeSyntheticExcursion()
  out <- propagateUncertainty(exc$time, exc$vals,
                              changeFun = detectExcursionCore,
                              n.ens = 15,
                              event.yr = 500, event.window = 100,
                              ref.window = 200)
  expect_equal(out, ref$propagated)
})

test_that("testNullHypothesis matches pre-migration reference", {
  exc <- makeSyntheticExcursion()
  nullOut <- testNullHypothesis(exc$time, exc$vals,
                                changeFun = detectExcursionCore,
                                n.ens = 5, mc.ens = 5,
                                surrogate.method = "isospectral",
                                event.yr = 500, event.window = 100,
                                ref.window = 200)
  expect_equal(purrr::map(nullOut, "eventDetected"), ref$nullEventDetected)
})

test_that("detectExcursion end-to-end matches pre-migration reference", {
  exc <- makeSyntheticExcursion()
  out <- detectExcursion(time = exc$time, vals = exc$vals,
                         time.units = "yr BP", vals.units = "unitless",
                         time.variable.name = "age",
                         vals.variable.name = "synthetic",
                         dataset.name = "syntheticExcursion",
                         n.ens = 15, null.hypothesis.n = 10,
                         event.yr = 500, event.window = 100,
                         ref.window = 200)
  expect_equal(out, ref$excursion)
  expect_s3_class(out, "excursion")
})

test_that("detectShift end-to-end matches pre-migration reference", {
  shf <- makeSyntheticShift()
  out <- detectShift(time = shf$time, vals = shf$vals,
                     time.units = "yr BP", vals.units = "unitless",
                     time.variable.name = "age",
                     vals.variable.name = "synthetic",
                     dataset.name = "syntheticShift",
                     n.ens = 15, null.hypothesis.n = 10,
                     summary.bin.step = 200,
                     minimum.segment.length = 100)
  expect_equal(out, ref$shift)
  expect_s3_class(out, "shift")
})
