#' Print EWS Detection Results
#'
#' @export
#' @param x \[`ews`]\cr An EWS detection result from [detect_warnings()].
#' @param ... Additional arguments passed to the generic print method.
#' @return `x` (invisibly).
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' ews <- detect_warnings(ts_data)
#' print(ews)
#'
print.ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print Regime Detection Results
#'
#' @export
#' @param x \[`regimes`]\cr A regime detection result from [detect_regimes()].
#' @param ... Additional arguments passed to the generic print method.
#' @return `x` (invisibly).
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' regimes <- detect_regimes(
#'   data = ts_data,
#'   method = "threshold",
#'   sensitivity = "medium"
#' )
#' print(regimes)
#'
print.regimes <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print Discovered Patterns
#'
#' @export
#' @param x \[`patterns`]\cr A pattern discovery result from
#'   [discover_patterns()].
#' @param ... Additional arguments passed to the generic print method.
#' @return `x` (invisibly).
#' @examples
#' ngrams <- discover_patterns(engagement, type = "ngram")
#' print(ngrams)
#'
print.patterns <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}
