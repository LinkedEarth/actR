# actR 0.4.0

* `prepareInput()`, and so `detectShift()` and `detectExcursion()`, use the time metadata that lipdR 0.7 adds in `extractTs()`. The time axis is an age or year ensemble when there is one, then lipdR's `time` (the primary time column), then `age`, then `year`. `timeVariableName` is now filled in from the time direction, so lipdR timeseries no longer need one. Requires lipdR >= 0.7.3.
* With a timeseries supplied, `time.units` now converts its time axis to that scale (e.g. `time.units = "yr BP"`), using `lipdR::getTsTime()`. `time.variable.name` can name the column to use (e.g. `"year"` or `"ageEnsemble"`).
* `detectMultipleExcursions()` gains `time.units` to put every record on one scale, so `event.yr` means the same time for all of them. Without it, records in different units now give a warning.
* `prepareInput()` and `detectMultipleExcursions()` accept a TS list (`lipdR::extractTs()`) and tibbles that have lost their class, as well as timeseries tibbles.
* `prepareInput()` no longer errors when a timeseries has no `paleoData_units` field; it fills in "?" as it does for missing values.
* Plots now reverse ka (and Ma, b2k) time axes, not just yr BP and kyr.

# actR 0.3.0

* actR now builds on ens and lipdViz instead of the geoChronR monolith, so age modeling software is no longer installed with it. Its uncertainty engine now comes from ens.
* Fixed a syntax error in `inst/CITATION`.
* Requires lipdR >= 0.7.0.

# actR 0.2.2
* Update detectShift() vignette following implementing `gaussianize = TRUE` as the default method.
* Fix bug in `prepareTimeseries()` Thanks again Georgy: [issue 32](https://github.com/LinkedEarth/actR/issues/32)

# actR 0.2.1

* Fixed bug in `plotExcursionSliding()`. Thanks Georgy [issue 31](https://github.com/LinkedEarth/actR/issues/31)

# actR 0.2.0

* Integrated some new features and plotting - most notably multisite and spatial significance testing, as well as mapping.

# actR 0.1.10

* Attempted to fix the bug that resulted from the rEDM v1.15.0 name change of `make_surrogate_data` to `SurrogateData`

# actR 0.1.9

* Excursion methods now export p-values for above, below, either or both options in all cases, to avoid unnecessary rerunning

# actR 0.1.8

* Add `detectExcursionSlidingWindow()`

# actR 0.1.7

* Introduce `detectMultipleExcursions()` to run the excursion detector across a suite of datasets
* Makes `detectExcursion()` return a data.frame with NAs rather than errors to better support `detectMultipleExcursions()`

# actR 0.1.6

* Used purrr >v1.0.0 to enable progress bars!
* fixed bug with exc.type and plotting
* minor fixes

# actR 0.1.5

* merged in calculate deltas - introducing the `summarizeChangepointMeanChanges()`
* More fixes for R >=4.2.0 sensitivity to multiple conditionals.

# actR 0.1.4

* Fix timeUnits bug in `prepareInput()` 
* Update `time.range` option for R >=4.2.0 sensitivity to multiple conditionals.

# actR 0.1.3

* added a `time.range` option to `prepareInput()` and `detectShift()` to optionally restrict the range of the analysis. Thanks to [@AndreaLemur](https://github.com/AndreaLemur) for the suggestion. 

# actR 0.1.2

* add `gaussianize` options to `detectShiftCore()`, with a default of TRUE, following [@AndreaLemur](https://github.com/AndreaLemur) recognition that that is a requirement for many of the methods in the `changepoint` package. More details on [changepoint workshop here](https://www.youtube.com/watch?v=UfGrLJ7S3sc) 

# actR 0.1.1

* added option to pass plot options to geoChronR plotting to plot.shift()
* Added better error handling with excursions

# actR 0.1.0

* This marks the first release of the package! Check here for updates
* Added a `NEWS.md` file to track changes to the package.
