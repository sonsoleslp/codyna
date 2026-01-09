prepare_sequence_data <- function(x, alphabet, cols) {
  check_missing(x)
  stopifnot_(
    is.data.frame(x) || is.matrix(x) || inherits(x, "stslist"),
    "Argument {.arg data} must be a {.cls data.frame}, a {.cls matrix} or
     an {.cls stslist} object."
  )
  p <- ncol(x)
  cols <- cols %m% seq_len(p)
  cols <- get_cols(rlang::enquo(cols), x)
  if (inherits(x, "stslist")) {
    alphabet <- attr(x, "alphabet")
    x <- as.data.frame(x)
  } else if (is.data.frame(x)) {
    if (missing(alphabet)) {
      vals <- sort(unique(unlist(x[, cols])))
      alphabet <- vals[!is.na(vals) & nchar(vals) > 0]
    }
    x[, cols] <- as.data.frame(
      lapply(x[, cols], function(y) factor(y, levels = alphabet))
    )
  } else if (inherits(x, "tna")) {
    if (!is.null(x$data)) {
      return(x$data)
    }
  }
  x <- as.matrix(
    as.data.frame(
      lapply(
        x[, cols],
        function(y) {
          as.integer(replace(y, which(!y %in% alphabet), NA))
        }
      )
    )
  )
  list(
    sequences = x,
    alphabet = alphabet
  )
}

prepare_timeseries_data <- function(x) {
  values <- as.numeric(x)
  time <- seq_along(values)
  if (stats::is.ts(x)) {
    tsp <- attr(x, "tsp")
    time <- seq(tsp[1], tsp[2], tsp[3])
  }
  list(values = values, time = time)
}
