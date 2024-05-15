# actR 0.2.0

* Integrated some new features and plotting - most notably multi site and spatial significance testing, as well as mapping.

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
