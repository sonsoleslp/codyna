#' Calculate Dynamic Complexity Measures for Time-Series Data
#'
#' Computes dynamic complexity and other rolling window measures for univariate
#' time series data.
#'
#' @export
#' @param data \[`ts`, `numeric()`]\cr
#'   Univariate time series data.
#' @param measures \[`character()`]\cr
#'   A vector of measures to calculate. See 'Details' for more information
#'   on the available measures.
#' @param window \[`integer(1)`: `7L`]\cr
#'   A positive `integer` specifying the rolling window size.
#'   Must be at least `2`.
#' @param align \[`character(1)`: `"center"`]\cr
#'   Alignment of the window. The available options are: `"center"`, `"right"`,
#'   and `"left"`. The calculated measure is assigned to the center, rightmost,
#'   or leftmost point of the window, respectively.
#' @return A `tibble` with the time index, the original time-series data,
#'   and the calculated measures.
#' @details The following measures can be calculated:
#'
#'   * `"complexity"`: Product of fluctuation and distribution measures.
#'   * `"fluctuation"`: Root mean square of successive differences.
#'   * `"distribution"`: Deviation from uniform distribution.
#'   * `"autocorrelation"`: Lag-1 autocorrelation coefficient.
#'   * `"max"`: Rolling maximum.
#'   * `"min"`: Rolling minimum.
#'   * `"variance"`: Rolling variance.
#'
#'   The option `"all"` computes all of the above.
#'
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#'
#' # Single measure
#' comp_single <- complexity(ts_data, measures = "complexity")
#'
#' # Multiple measures
#' comp_multi <- complexity(ts_data, measures = c("complexity", "variance"))
#'
complexity <- function(data, measures = "complexity",
                       window = 7L, align = "center") {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  valid_measures <- c(names(complexity_funs), "complexity")
  measures <- check_match(
    measures,
    c(valid_measures, "all"),
    several.ok = TRUE
  )
  measures <- ifelse_("all" %in% measures, valid_measures, measures)
  measures <- ifelse_(
    "complexity" %in% measures,
    c(setdiff(measures, "complexity"), "complexity"),
    measures
  )
  align <- check_match(align, c("left", "right", "center"))
  n <- length(values)
  check_range(window, type = "integer", min = 2L, max = n)
  scale <- range(values, na.rm = TRUE)
  out <- data.frame(values = values, time = time)
  for (measure in measures) {
    if (measure == "complexity") {
      fluctuation <- out$fluctuation %||% roll(
        fun = complexity_funs$fluctuation,
        values = values,
        window = window,
        align = align,
        scale = scale
      )
      distribution <- out$distribution %||% roll(
        fun = complexity_funs$distribution,
        values = values,
        window = window,
        align = align,
        scale = scale
      )
      out$complexity <- fluctuation * distribution
    } else {
      out[[measure]] <- roll(
        fun = complexity_funs[[measure]],
        values = values,
        window = window,
        align = align,
        scale = scale
      )
    }
  }
  tibble::as_tibble(out)
}

# Complexity methods ------------------------------------------------------

complexity_funs <- list(
  fluctuation = function(x, scale) {
    f_max <- scale[2L] - scale[1L]
    f_obs <- rmsqd(x)
    max(0, min(1, f_obs / f_max))
  },
  distribution = function(x, scale) {
    x <- x[!is.na(x)]
    n <- length(x)
    uniform <- seq(from = scale[1L], to = scale[2L], length.out = n)
    empirical <- sort(x)
    uni_diff <- diff(uniform)
    emp_diff <- diff(empirical)
    deviation <- uni_diff - emp_diff
    dev_h <- deviation * (sign(deviation) + 1) / 2
    div_diff <- dev_h / uni_diff
    div_diff[is.infinite(div_diff)] <- NA
    D <- 1 - mean(div_diff, na.rm = TRUE)
    max(0, min(1, D))
  },
  autocorrelation = function(x, ...) {
    x <- x[!is.na(x)]
    stats::acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2L]
  },
  min = function(x, ...) min(x, na.rm = TRUE),
  max = function(x, ...) max(x, na.rm = TRUE),
  variance = function(x, ...) stats::var(x, na.rm = TRUE)
)
