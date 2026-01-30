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
#' @param id \[`expression`]\cr A tidy selection of column names that
#'   uniquely identify each observation (optional).
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
convert <- function(data, cols = tidyselect::everything(), id,
                    format = "frequency") {
  check_missing(data)
  cols <- get_cols(rlang::enquo(cols), data)
  id <- get_cols(rlang::enquo(id), data)
  data <- prepare_sequence_data(data, cols, id)
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

prepare_sequence_data <- function(data, cols, id, group) {
  check_missing(data)
  stopifnot_(
    is.data.frame(data) ||
      inherits(data, "matrix") ||
      inherits(data, "stslist") ||
      inherits(data, "tna") ||
      inherits(data, "group_tna"),
    "Argument {.arg data} must be a {.cls data.frame}, a {.cls matrix},
     an {.cls stslist} object, a {.cls tna} object, or a {.cls group_tna}
    object."
  )
  parsed <- data_parsers[[class(data)[1L]]](data, cols)
  data <- parsed$data
  alphabet <- parsed$alphabet
  id <- id %m% NULL
  group <- group %m% NULL
  group <- ifelse_(is.null(parsed$group), group, parsed$group)
  n_group <- length(group)
  stopifnot_(
    is.null(group) || n_group == nrow(data) || n_group == 1L,
    "Argument {.arg group} must be either {.val last_obs}, a column name of
     {.arg data} or a {.cls vector} with the same length as the
     number of rows/sequences of {.arg data}."
  )
  if (n_group == 1L && group != "last_obs") {
    stopifnot_(
      group %in% names(data),
      "The column {.val {group}} was not found in the data."
    )
    tmp <- data[[group]]
    data[[group]] <- NULL
    group <- tmp
    cols <- setdiff(cols, group)
  }
  if (!is.null(id)) {
    cols <- setdiff(cols, id)
    id <- data[, id, drop = FALSE] |>
      as.data.frame() |>
      interaction(drop = TRUE)
  }
  if (is.data.frame(data)) {
    data <- data[, cols] |>
      lapply(
        function(y) {
          as.integer(replace(y, which(!y %in% alphabet), NA))
        }
      ) |>
      as.data.frame() |>
      as.matrix()
  }
  if (n_group == 1L && group == "last_obs") {
    nas <- is.na(data)
    last_obs <- max.col(!nas, ties.method = "last")
    group <- alphabet[data[, last_obs]]
    data[, last_obs] <- NA
    alphabet <- setdiff(alphabet, unique(group))
    stopifnot_(
      all(!x %in% alphabet),
      "Group identifiers must not be states of the sequence data."
    )
  }
  list(
    sequences = data,
    alphabet = alphabet,
    id = id,
    group = group
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


# Data parsers ------------------------------------------------------------


parse_tna <- function(x, ...) {
  stopifnot_(
    !is.null(x$data),
    "Argument {.arg data} is a {.cls tna} object with no data."
  )
  alphabet <- attr(x$data, "alphabet")
  out <- c(x$data)
  dim(out) <- dim(x$data)
  colnames(out) <- colnames(x$data)
  list(data = out, alphabet = alphabet)
}

parse_matrix <- function(x, cols) {
  stopifnot_(
    !is.null(colnames(x)),
    "Argument {.arg data} must have column names when a {.cls matrix} is
     provided."
  )
  parse_data.frame(as.data.frame(x), cols)
}

parse_group_tna <- function(x, ...) {
  cols <- attr(x, "cols")
  group <- attr(x, "groups")
  alphabet <- attr(x[[1L]]$data, "alphabet")
  data <- dplyr::bind_rows(
    lapply(x, function(y) as.data.frame(y$data))
  )
  group <- attr(x, "levels")[unlist(groups)]
  list(data = data, alphabet = alphabet, group = group)
}

parse_stslist <- function(x, ...) {
  list(data = as.data.frame(x), alphabet = attr(x, "alphabet"))
}

parse_data.frame <- function(x, cols) {
  vals <- as.character(sort(unique(unlist(x[, cols]))))
  alphabet <- vals[!is.na(vals) & nchar(vals) > 0]
  x[, cols] <- lapply(x[, cols], function(y) factor(y, levels = alphabet)) |>
    as.data.frame()
  list(data = x, alphabet = alphabet)
}

data_parsers <- list(
  tna = parse_tna,
  group_tna = parse_group_tna,
  matrix = parse_matrix,
  data.frame = parse_data.frame,
  stslist = parse_stslist
)
