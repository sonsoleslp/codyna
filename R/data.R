#' Convert Sequence Data to Various Formats
#'
#' Converts wide format sequence data into useful formats for analysis,
#' such as frequency table, one-Hot encoding, or edge list (graph format).
#'
#' @export
#' @param data \[`data.frame`]\cr
#'   Sequence data in wide format (rows are sequences, columns are time points).
#'   The input should be coercible to a `data.frame` object.
#' @param cols \[`expression`: `tidyselect::everything()`]\cr
#'   A tidy selection of columns that should be considered as sequence data.
#'   By default, all columns are used.
#' @param format \[`character(1)`: `"frequency"`]\cr
#'   The data format to convert into:
#'
#'   * `"frequency"`: Counts of each state per sequence.
#'   * `"onehot"`: Presence/absence (1/0) of each state per sequence.
#'   * `"edgelist"`: (state, next state) pairs.
#'   * `"reverse"`: Same as `"edgelist"` but in the reverse direction, i.e.,
#'     (state, previous state) pairs.
#'
#' @return A `tibble` structured according to the requested format.
#' @examples
#' convert(engagement, format = "frequency")
#' convert(engagement, format = "onehot")
#' convert(engagement, format = "edgelist")
#' convert(engagement, format = "reverse")
#'
convert <- function(data, cols = tidyselect::everything(),
                    format = "frequency") {
  check_missing(data)
  data <- extract_data(data)
  cols <- get_cols(rlang::enquo(cols), data)
  data <- prepare_sequence_data(data, cols)
  format <- check_match(
    format, c("frequency", "onehot", "edgelist", "reverse")
  )
  alphabet <- data$alphabet
  sequences <- as.data.frame(data$sequences)
  long <- sequences |>
    dplyr::mutate(.id = dplyr::row_number()) |>
    tidyr::pivot_longer(
      !(!!rlang::sym(".id")),
      values_to = "state",
      values_drop_na = TRUE
    ) |>
    dplyr::mutate(
      state = alphabet[!!rlang::sym("state")]
    )
  out <- NULL
  if (format == "frequency") {
    out <- long |>
      dplyr::count(!!rlang::sym(".id"), !!rlang::sym("state")) |>
      tidyr::pivot_wider(
        names_from = !!rlang::sym("state"),
        values_from = !!rlang::sym("n"),
        values_fill = 0
      )
  } else if (format == "onehot") {
    out <- long |>
      dplyr::distinct(!!rlang::sym(".id"), !!rlang::sym("state")) |>
      dplyr::mutate(present = 1) |>
      tidyr::pivot_wider(
        names_from = !!rlang::sym("state"),
        values_from = !!rlang::sym("present"),
        values_fill = 0
      )
  } else if (format == "edgelist") {
    out <- long |>
      dplyr::group_by(!!rlang::sym(".id")) |>
      dplyr::mutate(to = dplyr::lead(!!rlang::sym("state"))) |>
      dplyr::filter(!is.na(!!rlang::sym("to"))) |>
      dplyr::select(
        !!rlang::sym(".id"),
        from = !!rlang::sym("state"),
        !!rlang::sym("to")
      ) |>
      dplyr::ungroup()
  } else if (format == "reverse") {
    out <- long |>
      dplyr::group_by(!!rlang::sym(".id")) |>
      dplyr::mutate(previous = dplyr::lag(!!rlang::sym("state"))) |>
      dplyr::filter(!is.na(!!rlang::sym("previous"))) |>
      dplyr::select(
        !!rlang::sym(".id"),
        !!rlang::sym("state"),
        !!rlang::sym("previous")
      ) |>
      dplyr::ungroup()
  }
  out
}

prepare_sequence_data <- function(x, cols) {
  alphabet <- attr(x, "alphabet")
  if (is.null(alphabet)) {
    vals <- as.character(sort(unique(unlist(x[, cols]))))
    alphabet <- vals[!is.na(vals) & nchar(vals) > 0]
  }
  x <- x[, cols] |>
    lapply(
      function(y) {
        y <- factor(y, levels = alphabet)
        as.integer(replace(y, which(!y %in% alphabet), NA))
      }
    ) |>
    as.data.frame() |>
    as.matrix()
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
    time <- seq(tsp[1L], tsp[2L], tsp[3L])
  }
  list(values = values, time = time)
}

extract_data <- function(x, group = FALSE, meta = FALSE) {
  if (is.matrix(x)) {
    stopifnot_(
      !is.null(colnames(x)),
      "Argument {.arg data} must have column names when a {.cls matrix} is
       provided."
    )
    return(as.data.frame(x))
  }
  if (inherits(x, "tna")) {
    stopifnot_(
      !is.null(x$data),
      "Argument {.arg data} is a {.cls tna} object with no data."
    )
    alphabet <- attr(x$data, "alphabet")
    out <- alphabet[c(x$data)]
    dim(out) <- dim(x$data)
    colnames(out) <- colnames(x$data)
    out <- as.data.frame(out)
    attr(out, "alphabet") <- attr(x$data, "alphabet")
    return(out)
  }
  if (inherits(x, "group_tna")) {
    alphabet <- attr(x[[1L]]$data, "alphabet")
    data <- do.call(base::rbind, lapply(x, "[[", "data"))
    out <- alphabet[c(data)]
    dim(out) <- dim(data)
    colnames(out) <- colnames(data)
    out <- as.data.frame(out)
    attr(out, "alphabet") <- alphabet
    if (group) {
      group <- attr(x, "groups")
      out$.group <- attr(x, "levels")[unlist(group)]
      attr(out, "group") <- ".group"
    }
    return(out)
  }
  if (inherits(x, "tna_data")) {
    out <- x$sequence_data
    attr(out, "cols") <- names(out)
    if (meta) {
      out <- cbind(x$meta_data, out)
    }
    return(out)
  }
  x
}

extract_outcome <- function(x, outcome) {
  if (missing(outcome)) {
    return(list(last = FALSE, outcome = NULL, var = NULL))
  }
  n_out <- length(outcome)
  stopifnot_(
    n_out == nrow(x) || n_out == 1L,
    "Argument {.arg outcome} must be either {.val last_obs}, a column name of
     {.arg data} or a {.cls vector} with the same length as the
     number of rows of {.arg data}."
  )
  if (n_out == 1L) {
    if (outcome == "last_obs") {
      return(list(last = TRUE, outcome = NULL, var = NULL))
    }
    outcome <- as.character(outcome)
    stopifnot_(
      outcome %in% names(x),
      "The column {.val {outcome}} must exist in the data."
    )
    return(list(last = FALSE, outcome = x[[outcome]]), var = outcome)
  }
  list(last = FALSE, outcome = outcome, var = NULL)
}

extract_last <- function(x, alphabet) {
  n <- nrow(x)
  nas <- is.na(x)
  last_obs <- max.col(!nas, ties.method = "last")
  idx <- cbind(seq_len(n), last_obs)
  last <- x[idx]
  last_vals <- unique(last)
  group <- alphabet[last]
  x[idx] <- NA
  x[x %in% last_vals] <- NA
  groups <- unique(group)
  vals <- seq_along(alphabet)
  alphabet <- setdiff(alphabet, groups)
  vals[last_vals] <- NA
  vals[-last_vals] <- seq_along(alphabet)
  d <- dim(x)
  x <- vals[x]
  dim(x) <- d
  list(sequences = x, alphabet = alphabet, group = group)
}
