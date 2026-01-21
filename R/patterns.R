#' Discover Sequence Patterns
#'
#' Discovering various types of patterns in sequence data.
#' Provides n-gram extraction, gapped pattern discovery, analysis of repeated
#' patterns and targeted pattern search.
#'
#' @export
#' @param data \[`data.frame`, `matrix`, `stslist`]\cr
#'   Sequence data in wide format (rows are sequences, columns are time points).
#' @param type \[`character(1)`]\cr Type of pattern analysis:
#'
#'   * `"ngram"`: Extract contiguous n-grams.
#'   * `"gapped"`: Discover patterns with gaps/wildcards.
#'   * `"repeated"`: Detect repeated occurrences of the same state.
#'
#' @param pattern \[`character(1)`]\cr Specific pattern to search for as
#'   a character string (e.g., `"A->*->B"`). If provided, `type` is ignored.
#'   Supports wildcards: `*` (single) and `**` (multi-wildcard).
#' @param len \[`integer()`]\cr Pattern lengths to consider for
#'   n-grams and repeated patterns (default: `2:5`).
#' @param gap \[`integer()`]\cr Gap sizes to consider for
#'   gapped patterns (default: `1:3`).
#' @param min_support \[`integer(1)`]\cr Minimum support threshold, i.e., the
#'   proportion of sequences that must contain a specific pattern for it to
#'   be included (default: `0.01`).
#' @param min_count \[`integer(1)`]\cr Minimum count threshold, i.e., the
#'   numbers of times a pattern must occur across all sequences for it to
#'   be included (default: `2`).
#' @param start \[`character(1)`]\cr Filter patterns starting with these states.
#' @param end \[`character(1)`]\cr Filter patterns ending with these states.
#' @param contains \[`character(1)`]\cr Filter patterns containing these states.
#' @return A `tibble` containing the discover patterns, counts, proportions,
#'   support and lift.
#' @examples
#' # N-grams
#' ngrams <- discover_patterns(engagement, type = "ngram")
#'
#' # Gapped patterns
#' gapped <- discover_patterns(engagement, type = "gapped")
#'
#' # Repeated patterns
#' repeated <- discover_patterns(engagement, type = "repeated")
#'
#' # Custom pattern with a wildcard state
#' custom <- discover_patterns(engagement, pattern = "Active->*")
#'
discover_patterns <- function(data, type = "ngram", pattern, len = 2:5,
                              gap = 1:3, min_support = 0.01, min_count = 2,
                              start, end, contains) {
  check_missing(data)
  data <- prepare_sequence_data(data)
  sequences <- data$sequences
  alphabet <- data$alphabet
  m <- ncol(sequences)
  type <- check_match(type, c("ngram", "gapped", "repeated"))
  check_range(len, scalar = FALSE, type = "integer", min = 2L, max = m)
  check_range(gap, scalar = FALSE, type = "integer", min = 1L, max = m - 2L)
  check_range(min_support)
  check_values(min_count)
  if (!missing(pattern)) {
    check_string(pattern)
    patterns <- search_pattern(sequences, alphabet, pattern)
  } else {
    patterns <- switch(type,
      `ngram` = extract_ngrams(sequences, alphabet, len),
      `gapped` = extract_gapped(sequences, alphabet, gap),
      `repeated` = extract_repeated(sequences, alphabet, len)
    )
  }
  support <- state_support(sequences, alphabet)
  process_patterns(patterns) |>
    filter_patterns(min_support, min_count, start, end, contains) |>
    pattern_proportions() |>
    pattern_lift(support = support)
}

process_patterns <- function(x) {
  out <- vector(mode = "list", length = length(x))
  for (i in seq_along(out)) {
    pat_mat <- x[[i]]$patterns
    pat_len <- x[[i]]$length
    n <- nrow(pat_mat)
    pat_vec <- c(pat_mat[nzchar(pat_mat)])
    if (length(pat_vec) == 0) {
      next
    }
    pat_tab <- table(pat_vec)
    pat_u <- names(pat_tab)
    pat_con <- integer(length(pat_u))
    names(pat_con) <- pat_u
    for (j in seq_len(n)) {
      pats <- pat_mat[j, ]
      pats <- pats[nzchar(pats)]
      pat_con[pats] <- pat_con[pats] + 1L
    }
    out[[i]] <- data.frame(
      pattern = pat_u,
      length = pat_len,
      count = c(pat_tab),
      contained_in = pat_con,
      support = pat_con / n
    )
  }
  proto <- data.frame(
    pattern = character(0L),
    length = integer(0L),
    count = integer(0L),
    contained_in = integer(0L),
    support = numeric(0L)
  )
  out <- dplyr::bind_rows(out, proto) |>
    dplyr::arrange(dplyr::desc(!!rlang::sym("count"))) |>
    tibble::as_tibble()
}

extract_ngrams <- function(sequences, alphabet, len) {
  n <- nrow(sequences)
  m <- ncol(sequences)
  k <- length(len)
  mis <- is.na(sequences)
  ngrams <- vector(mode = "list", length = k)
  for (l in seq_len(k)) {
    j <- len[l]
    tmp <- matrix("", n, m - j + 1L)
    for (i in seq_len(m - j + 1L)) {
      idx <- i:(i + j - 1L)
      pattern <- alphabet[sequences[, idx]]
      dim(pattern) <- c(n, j)
      pattern_mis <- mis[, idx]
      valid <- .rowSums(pattern_mis, m = n, n = j) == 0L
      if (any(valid)) {
        subseq <- do.call(
          paste,
          c(
            as.data.frame(pattern[valid, , drop = FALSE]),
            sep = "->"
          )
        )
        tmp[valid, i] <- subseq
      }
    }
    ngrams[[l]] <- list(patterns = tmp, length = j)
  }
  ngrams
}

extract_gapped <- function(sequences, alphabet, gap) {
  n <- nrow(sequences)
  m <- ncol(sequences)
  k <- length(gap)
  mis <- is.na(sequences)
  gapped <- vector(mode = "list", length = k)
  for (l in seq_len(k)) {
    j <- gap[l]
    tmp <- matrix("", n, m - j)
    for (i in seq_len(m - j - 1L)) {
      idx <- c(i, i + j + 1L)
      pattern <- alphabet[sequences[, idx]]
      dim(pattern) <- c(n, 2L)
      pattern_mis <- mis[, idx]
      valid <- .rowSums(pattern_mis, m = n, n = 2L) == 0L
      if (any(valid)) {
        subseq <- do.call(
          paste,
          c(
            as.data.frame(pattern[valid, , drop = FALSE]),
            sep = paste0("->", paste0(rep("*", j), collapse = ""), "->")
          )
        )
        tmp[valid, i] <- subseq
      }
    }
    gapped[[l]] <- list(patterns = tmp, length = j + 2L)
  }
  gapped
}

extract_repeated <- function(sequences, alphabet, len) {
  n <- nrow(sequences)
  m <- ncol(sequences)
  k <- length(len)
  mis <- is.na(sequences)
  repeated <- vector(mode = "list", length = k)
  for (l in seq_len(k)) {
    j <- len[l]
    tmp <- matrix("", n, m - j + 1L)
    for (i in seq_len(m - j + 1L)) {
      idx <- i:(i + j - 1L)
      pattern <- alphabet[sequences[, idx]]
      dim(pattern) <- c(n, j)
      pattern_mis <- mis[, idx]
      valid <- (.rowSums(pattern_mis, m = n, n = j) == 0L) &
        apply(sequences[, idx], 1, function(x) n_unique(x) == 1L)
      if (any(valid)) {
        subseq <- do.call(
          paste,
          c(
            as.data.frame(pattern[valid, , drop = FALSE]),
            sep = "->"
          )
        )
        tmp[valid, i] <- subseq
      }
    }
    repeated[[l]] <- list(patterns = tmp, length = j)
  }
  repeated
}

search_pattern <- function(sequences, alphabet, pattern) {
  n <- nrow(sequences)
  m <- ncol(sequences)
  mis <- is.na(sequences)
  states <- strsplit(pattern, split = "->")[[1]]
  wildcards <- grepl("^\\*+$", states, perl = TRUE)
  j <- length(states)
  pos <- seq_len(j)
  if (any(wildcards)) {
    idx_wild <- which(wildcards)
    idx_states <- which(!wildcards)
    s <- length(idx_states)
    pos <- integer(s)
    for (i in seq_len(s)) {
      k <- idx_states[i]
      pos[i] <- i + sum(nchar(states[idx_wild[idx_wild < k]]))
    }
    j <- s + sum(nchar(states[idx_wild]))
    states <- states[idx_states]
  }
  if (j > m) {
    return(list(list(patterns = matrix("", nrow = 0, ncol = j), length = j)))
  }
  discovered <- matrix("", n, m - j + 1L)
  for (i in seq_len(m - j + 1L)) {
    idx <- i:(i + j - 1L)
    pattern <- alphabet[sequences[, idx]]
    dim(pattern) <- c(n, j)
    pattern_mis <- mis[, idx]
    valid <- (.rowSums(pattern_mis, m = n, n = j) == 0L) &
      apply(pattern, 1, function(x) all(x[pos] == states))
    if (any(valid)) {
      subseq <- do.call(
        paste,
        c(
          as.data.frame(pattern[valid, , drop = FALSE]),
          sep = "->"
        )
      )
      discovered[valid, i] <- subseq
    }
  }
  list(list(patterns = discovered, length = j))
}

filter_patterns <- function(patterns, min_support, min_count,
                            start, end, contains) {
  out <- patterns |>
    dplyr::filter(!!rlang::sym("support") >= min_support) |>
    dplyr::filter(!!rlang::sym("count") >= min_count)
  if (!missing(start)) {
    is_prefix <- function(x, prefix) {
      vapply(x, function(y) any(startsWith(y, prefix)), logical(1L))
    }
    out <- out |>
      dplyr::filter(is_prefix(!!rlang::sym("pattern"), start))
  }
  if (!missing(end)) {
    is_suffix <- function(x, suffix) {
      vapply(x, function(y) any(endsWith(y, suffix)), logical(1L))
    }
    out <- out |>
      dplyr::filter(is_suffix(!!rlang::sym("pattern"), end))
  }
  if (!missing(contains)) {
    pat <- paste0(contains, collapse = "|")
    out <- out |>
      dplyr::filter(grepl(pat, !!rlang::sym("pattern"), perl = TRUE))
  }
  out
}

state_support <- function(sequences, alphabet) {
  n <- nrow(sequences)
  m <- ncol(sequences)
  a <- length(alphabet)
  support <- integer(a)
  for (i in seq_len(a)) {
    support[i] <- sum(
      .rowSums(sequences == i, m = n, n = m, na.rm = TRUE) > 0
    )
  }
  names(support) <- alphabet
  support / n
}

pattern_proportions <- function(patterns) {
  patterns |>
    dplyr::group_by(!!rlang::sym("length")) |>
    dplyr::mutate(
      proportion = !!rlang::sym("count") / sum(!!rlang::sym("count"))
    ) |>
    dplyr::relocate(
      !!rlang::sym("proportion"),
      .after = !!rlang::sym("count")
    ) |>
    dplyr::ungroup()
}

pattern_lift <- function(patterns, support) {
  patterns$lift <- NA_real_
  pattern_states <- strsplit(patterns$pattern, split = "->")
  n <- nrow(patterns)
  denom <- numeric(n)
  for (i in seq_len(n)) {
    states <- pattern_states[[i]]
    states <- states[!grepl("\\*+", states)]
    denom[i] <- prod(support[states])
  }
  patterns$lift <- patterns$support / denom
  patterns
}
