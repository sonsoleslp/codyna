#' Convert Sequence Data to Various Formats
#'
#' Converts wide format sequence data into useful formats for analysis,
#' such as frequency table, one-Hot encoding, or edge list (graph format).
#'
#' @export
#' @param data \[`data.frame`, `matrix`, `stslist`]\cr
#'   Sequence data in wide format (rows are sequences, columns are time points).
#' @param cols \[`expression`]\cr A tidy selection of columns that should
#'   be considered as sequence data. By default, all columns are used.
#' @param format \[`character(1)`]\cr The format to convert into:
#'
#'   * `"frequency"`: Counts of each state per sequence.
#'   * `"onehot"`: Presence/absence (0/1) of each state per sequence.
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
convert <- function(data, cols, format = "frequency") {
  check_missing(data)
  data <- prepare_sequence_data(data, cols)
  format <- check_match(format, c("frequency", "onehot", "edgelist", "reverse"))
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

prepare_sequence_data <- function(x, alphabet, cols, group) {
  check_missing(x)
  stopifnot_(
    is.data.frame(x) ||
      inherits(x, "stslist") ||
      inherits(x, "tna"),
    "Argument {.arg data} must be a {.cls data.frame}, a {.cls matrix},
     an {.cls stslist} object, or a {.cls tna} object."
  )
  if (inherits(x, "tna")) {
    if (!is.null(x$data)) {
      sequences <- x$data
      alphabet <- attr(x$data, "alphabet")
      attr(sequences, "labels") <- NULL
      attr(sequences, "colors") <- NULL
      attr(sequences, "alphabet") <- NULL
      class(sequences) <- "matrix"
      return(list(sequences = sequences, alphabet = alphabet))
    }
    stop_("Argument {.arg data} is a {.cls tna} object with no data.")
  }
  group <- group %m% NULL
  n_group <- length(group)
  stopifnot_(
    is.null(group) || n_group == nrow(x) || n_group == 1L,
    "Argument {.arg group} must be either {.val last_obs}, a column name of
     {.arg data} or a {.cls vector} with the same length as the
     number of rows/sequences of {.arg data}."
  )
  if (n_group == 1L && group != "last_obs"){
    check_cols(group, names(x))
    tmp <- x[[group]]
    x[[group]] <- NULL
    group <- tmp
  }
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
  if (n_group == 1L && group == "last_obs") {
    nas <- is.na(x)
    last_obs <- max.col(!nas, ties.method = "last")
    group <- alphabet[x[, last_obs]]
    x[, last_obs] <- NA
    alphabet <- setdiff(alphabet, unique(group))
    stopifnot_(
      all(!x %in% alphabet),
      "Group identifiers must not be states of the sequence data."
    )
  }
  list(
    sequences = x,
    alphabet = alphabet,
    group = group
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
