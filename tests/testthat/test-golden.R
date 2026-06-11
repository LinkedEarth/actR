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

# Surrogate generation (rEDM::SurrogateData, used by testNullHypothesis) is
# seed-stable within a session but NOT reproducible across R/rEDM versions, so
# null-distribution and p-value columns cannot be compared bit-exactly between
# the capture environment and CI. For the surrogate-dependent results we compare
# the deterministic columns exactly and validity-check the stochastic ones.
# surrogate-derived columns whose VALUES are not portable across R/rEDM
# versions. null.hypothesis.n / unc.prop.n are deterministic counts and stay in
# the exact-comparison set.
isStochasticName <- function(nm) grepl("null_probability|^pvalue", nm, ignore.case = TRUE)
isProbabilityName <- function(nm) grepl("null_probability|^pvalue", nm, ignore.case = TRUE)

expect_equal_deterministic <- function(out, ref) {
  shared <- intersect(names(out), names(ref))
  det <- shared[!isStochasticName(shared)]
  for (nm in det) expect_equal(out[[nm]], ref[[nm]], info = nm)
}

expect_valid_probabilities <- function(out) {
  pcols <- names(out)[isProbabilityName(names(out))]
  for (nm in pcols) {
    v <- suppressWarnings(as.numeric(unlist(out[[nm]])))
    v <- v[is.finite(v)]
    if (length(v)) expect_true(all(v >= 0 & v <= 1), info = nm)
  }
}

test_that("testNullHypothesis returns well-formed per-surrogate detections", {
  exc <- makeSyntheticExcursion()
  nullOut <- testNullHypothesis(exc$time, exc$vals,
                                changeFun = detectExcursionCore,
                                n.ens = 5, mc.ens = 5,
                                surrogate.method = "isospectral",
                                event.yr = 500, event.window = 100,
                                ref.window = 200)
  # structure is portable even though surrogate values are not
  expect_length(nullOut, length(ref$nullEventDetected))
  detected <- purrr::map(nullOut, "eventDetected")
  expect_true(all(purrr::map_lgl(detected, ~ is.logical(.x) || is.numeric(.x))))
})

test_that("detectExcursion end-to-end matches reference (deterministic cols)", {
  exc <- makeSyntheticExcursion()
  out <- detectExcursion(time = exc$time, vals = exc$vals,
                         time.units = "yr BP", vals.units = "unitless",
                         time.variable.name = "age",
                         vals.variable.name = "synthetic",
                         dataset.name = "syntheticExcursion",
                         n.ens = 15, null.hypothesis.n = 10,
                         event.yr = 500, event.window = 100,
                         ref.window = 200)
  expect_s3_class(out, "excursion")
  expect_equal_deterministic(out, ref$excursion)
  expect_valid_probabilities(out)
})

test_that("detectShift end-to-end matches reference (deterministic cols)", {
  shf <- makeSyntheticShift()
  out <- detectShift(time = shf$time, vals = shf$vals,
                     time.units = "yr BP", vals.units = "unitless",
                     time.variable.name = "age",
                     vals.variable.name = "synthetic",
                     dataset.name = "syntheticShift",
                     n.ens = 15, null.hypothesis.n = 10,
                     summary.bin.step = 200,
                     minimum.segment.length = 100)
  expect_s3_class(out, "shift")
  expect_equal_deterministic(out, ref$shift)
  expect_valid_probabilities(out)
})
