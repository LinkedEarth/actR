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

# Some output columns embed values from third-party packages whose results are
# not reproducible across versions/platforms, so they cannot be compared
# bit-exactly between the capture environment (R 4.4 / macOS) and CI (R 4.6,
# Linux/Windows/macOS):
#   - null_probability* / pvalue* / cl*           : rEDM::SurrogateData
#   - parameters / method / penalty / pen.value /
#     ncpts.max / cpt.fun                          : changepoint::cpt internals
#                                                    (e.g. default method AMOC
#                                                    vs PELT, MBIC penalty value)
#   - *_hash (e.g. it_hash)                        : digest hashes of
#                                                    floating-point data; a 1-ulp
#                                                    difference across CPU
#                                                    architectures (macOS ARM vs
#                                                    x86) changes the whole hash
# Everything else -- time_start/end/mid, event_probability*, deltas, counts,
# minimum.segment.length, and metadata -- is portable. null.hypothesis.n /
# unc.prop.n are deterministic counts and remain in the exact-comparison set.
isStochasticName <- function(nm) {
  grepl(paste0("null_probability|^pvalue|^cl[0-9.]+$|^parameters$|",
               "^method$|^penalty$|^pen\\.value$|^ncpts|^cpt\\.fun$|hash$"),
        nm, ignore.case = TRUE)
}
isProbabilityName <- function(nm) grepl("null_probability|^pvalue", nm, ignore.case = TRUE)

expect_equal_deterministic <- function(out, ref) {
  shared <- intersect(names(out), names(ref))
  for (nm in shared) {
    if (isStochasticName(nm)) next
    if (is.list(out[[nm]]) && !is.data.frame(out[[nm]])) next # skip nested list-columns
    expect_equal(out[[nm]], ref[[nm]], info = nm)
  }
}

expect_valid_probabilities <- function(out) {
  pcols <- names(out)[isProbabilityName(names(out))]
  for (nm in pcols) {
    v <- suppressWarnings(as.numeric(unlist(out[[nm]])))
    v <- v[is.finite(v)]
    if (length(v)) expect_true(all(v >= 0 & v <= 1), info = nm)
  }
}

test_that("detectExcursionCore matches pre-migration reference", {
  exc <- makeSyntheticExcursion()
  out <- detectExcursionCore(exc$time, exc$vals,
                             event.yr = 500, event.window = 100,
                             ref.window = 200)
  expect_equal(out, ref$excursionCore)
  expect_true(out$eventDetected)
})

test_that("detectShiftCore matches pre-migration reference (deterministic cols)", {
  shf <- makeSyntheticShift()
  out <- detectShiftCore(shf$time, shf$vals, minimum.segment.length = 100)
  # the `parameters` string embeds changepoint::cpt internals (pen.value,
  # method) that vary across changepoint package versions; the numeric
  # detection columns are portable.
  expect_equal_deterministic(out, ref$shiftCore)
  expect_type(out$parameters, "character")
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
