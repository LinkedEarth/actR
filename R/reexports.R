# Generic ensemble-uncertainty machinery now lives in the ens package. actR
# pioneered this propagate -> null-test -> p-value pattern; it was promoted to
# ens so the whole package family (correlation, regression, compositing, etc.)
# can reuse it. These re-exports keep the functions available from actR so
# existing code and the internal detectors continue to work unchanged.

#' @importFrom ens surrogateDataFun
#' @export
ens::surrogateDataFun

#' @importFrom ens simulateAutoCorrelatedUncertainty
#' @export
ens::simulateAutoCorrelatedUncertainty

#' @importFrom ens propagateUncertainty
#' @export
ens::propagateUncertainty

#' @importFrom ens testNullHypothesis
#' @export
ens::testNullHypothesis

#' @importFrom ens kdePval
#' @export
ens::kdePval
