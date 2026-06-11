# Capture golden reference outputs from the CURRENT implementation of actR's
# core detection, uncertainty propagation, and null-hypothesis machinery,
# before rewiring dependencies to the ens/lipdViz stack. Run from package root.
devtools::load_all(".", quiet = TRUE)

ref <- list()

# deterministic synthetic series: excursion near t = 500
makeSyntheticExcursion <- function() {
  set.seed(42)
  time <- seq(0, 1000, by = 10)
  vals <- rnorm(length(time))
  vals[48:52] <- vals[48:52] + 4
  list(time = time, vals = vals)
}

# deterministic synthetic series: mean shift at t = 500
makeSyntheticShift <- function() {
  set.seed(43)
  time <- seq(0, 1000, by = 10)
  vals <- rnorm(length(time)) + ifelse(time > 500, 2, 0)
  list(time = time, vals = vals)
}

exc <- makeSyntheticExcursion()
shf <- makeSyntheticShift()

# 1. core detectors (fully deterministic)
ref$excursionCore <- detectExcursionCore(exc$time, exc$vals,
                                         event.yr = 500, event.window = 100,
                                         ref.window = 200)
ref$shiftCore <- detectShiftCore(shf$time, shf$vals,
                                 minimum.segment.length = 100)

# 2. uncertainty propagation (deterministic via default seed)
ref$propagated <- propagateUncertainty(exc$time, exc$vals,
                                       changeFun = detectExcursionCore,
                                       n.ens = 15,
                                       event.yr = 500, event.window = 100,
                                       ref.window = 200)

# 3. null hypothesis testing (deterministic via default seed)
nullOut <- testNullHypothesis(exc$time, exc$vals,
                              changeFun = detectExcursionCore,
                              n.ens = 5, mc.ens = 5,
                              surrogate.method = "isospectral",
                              event.yr = 500, event.window = 100,
                              ref.window = 200)
ref$nullEventDetected <- purrr::map(nullOut, "eventDetected")

# 4. end-to-end excursion detection (small ensembles for speed)
ref$excursion <- detectExcursion(time = exc$time, vals = exc$vals,
                                 time.units = "yr BP", vals.units = "unitless",
                                 time.variable.name = "age",
                                 vals.variable.name = "synthetic",
                                 dataset.name = "syntheticExcursion",
                                 n.ens = 15, null.hypothesis.n = 10,
                                 event.yr = 500, event.window = 100,
                                 ref.window = 200)

# 5. end-to-end shift detection
ref$shift <- detectShift(time = shf$time, vals = shf$vals,
                         time.units = "yr BP", vals.units = "unitless",
                         time.variable.name = "age",
                         vals.variable.name = "synthetic",
                         dataset.name = "syntheticShift",
                         n.ens = 15, null.hypothesis.n = 10,
                         summary.bin.step = 200,
                         minimum.segment.length = 100)

saveRDS(ref, "tests/testthat/reference_outputs.rds")
cat("Saved", length(ref), "reference objects\n")
