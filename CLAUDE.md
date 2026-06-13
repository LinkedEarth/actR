# CLAUDE.md — actR

`actR` (Abrupt Change Toolkit in R) detects, quantifies, and visualizes abrupt changes in
paleogeoscientific timeseries. As of June 2026 it is part of a six-package family and builds
directly on the **ens** engine (it no longer depends on the monolithic geoChronR).

Repo: **LinkedEarth/actR** (NOT nickmckay — `gh` API 404s on nickmckay/actR). Branch: `refactor`.

## Package family (dependency DAG: ens ← lipdViz ← geoChronR; actR, compositeR & fluxcapacitoR on top)

| Repo (`~/GitHub/...`) | GitHub | Branch | Role |
|---|---|---|---|
| ens | nickmckay/ens | main | Ensemble methods + UQ engine (used heavily here) |
| lipdViz | nickmckay/lipdViz | main | Plotting (ensemble ribbons) |
| geoChronR-chronOnly | nickmckay/geoChronR-chronOnly | main | Age modeling |
| **actR** (this repo) | LinkedEarth/actR | refactor | Abrupt-change detection |
| compositeR | nickmckay/compositeR | refactor | Record compositing |
| fluxcapacitoR | nickmckay/fluxcapacitoR | main | Flux-focused varve age modeling |

## Architecture

Each detector is a **`changeFun(time, vals, ...) -> one-row tibble`** (`detectExcursionCore`,
`detectShiftCore`). The generic engine that wraps them — `propagateUncertainty`,
`testNullHypothesis`, `kdePval`, `surrogateDataFun`, `simulateAutoCorrelatedUncertainty` —
**now lives in ens** and is re-exported here via `R/reexports.R`. To add a detector, write a
Core function; propagation + null-hypothesis testing come for free. Results are S3-classed
tibbles (`excursion`/`shift` + `*Core`).

## Gotchas

- Uses **markdown roxygen** (`Roxygen: list(markdown = TRUE)`); ens uses plain roxygen. Keep
  ens param text escape-clean — `@inheritParams ens::propagateUncertainty` re-emits it raw.
- **Golden tests** (`tests/testthat/test-golden.R`, ref via `tools/capture_reference.R`):
  many actR outputs embed values NOT reproducible across CI's R/pkg versions or macOS-ARM-vs-x86
  — rEDM surrogates (`null_probability*`/`pvalue*`/`cl*`), changepoint internals
  (`method`/`penalty`/`pen.value`/`ncpts.max`/`cpt.fun`/`parameters`), and **digest hashes of
  floats** (`*_hash`, e.g. `it_hash` — 1-ulp FP diff changes the whole hash). The test's
  `isStochasticName()` denylists all of these; only deterministic R-RNG numeric columns are
  compared exactly. Enumerate non-portable columns up front when adding golden tests.
- CI: R CMD check on Windows/Linux/macOS with vignettes. `gh run watch` is flaky — poll
  `gh run view <id> --json status` in an until-loop.

## Deferred TODO

- **vctrs / `dplyr_reconstruct` class-restore** so dplyr verbs (filter/mutate/select/slice)
  don't strip the `excursion`/`shift` S3 class. NOT a quick fix — needs vctrs proxy/restore +
  dplyr_reconstruct methods done carefully; a half-done version risks breaking tibble invariants.

## Dev

`devtools::load_all()` · `devtools::document()` · `devtools::test()` · `devtools::check()`.
Needs `ens` + `lipdViz` installed from source. Commit work when complete.
Co-author trailer: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
