# actR 0.1.6

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
