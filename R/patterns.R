# Discovery ---------------------------------------------------------------


#' Discover Sequence Patterns
#'
#' Discover various types of patterns in sequence data.
#' Provides n-gram extraction, gapped pattern discovery, analysis of repeated
#' patterns and targeted pattern search. supports comparison of pattern
#' presence between groups.
#'
#' @export
#' @inheritParams convert
#' @param group \[`character(1)`, `vector()`]\cr How should the sequences
#'   be grouped? The option `"last_obs"` assumes that the last non-missing
#'   observation of each sequence specifies the group membership. Alternatively,
#'   a column name of `data` or a `vector` with the same length as the number
#'   of rows. If not provided (default), the sequences will not not grouped.
#'   When provided, the presence/absence of each pattern is compared between
#'   the groups using a goodness-of-fit chi-square test.
#' @param type \[`character(1)`]\cr Type of pattern analysis:
#'
#'   * `"ngram"`: Extract contiguous n-grams.
#'   * `"gapped"`: Discover patterns with gaps/wildcards.
#'   * `"repeated"`: Detect repeated occurrences of the same state.
#'
#' @param pattern \[`character(1)`]\cr Specific pattern to search for as
#'   a character string (e.g., `"A->*->B"`). If provided, `type` is ignored.
#'   Supports wildcards: `*` (single) and `**` (multi-wildcard).
#' @param len \[`integer()`]\cr Pattern lengths to consider for n-grams and
#'   repeated patterns (default: `2:5`).
#' @param gap \[`integer()`]\cr Gap sizes to consider for gapped patterns
#'   (default: `1:3`).
#' @param min_support \[`integer(1)`]\cr Minimum support threshold, i.e., the
#'   proportion of sequences that must contain a specific pattern for it to
#'   be included (default: `0.01`).
#' @param min_freq \[`integer(1)`]\cr Minimum pattern frequency threshold, i.e.,
#'   the numbers of times a pattern must occur across all sequences for it to
#'   be included (default: `2`).
#' @param start \[`character(1)`]\cr Filter patterns starting with these states.
#' @param end \[`character(1)`]\cr Filter patterns ending with these states.
#' @param contains \[`character(1)`]\cr Filter patterns containing these states.
#' @return An object of class `patterns` which is a `tibble` with the
#'   following columns:
#'
#'   * `pattern`: the discovered patterns.
#'   * `length`: the length of the pattern.
#'   * `frequency`: the number of times the pattern occurs across all sequences.
#'   * `count`: the number of sequences the pattern occurs in.
#'   * `support`: the proportion of sequences that contain the pattern.
#'   * `lift`: the support divided by the product of the supports of the
#'     individual states of the pattern. For wildcards, the support is always 1.
#'
#' In addition, if `group` is provided, additional columns giving the counts
#' in each group, the chi-squared test statistic values, and p-values are
#' included. If `group` has two classes, the degree of separation is also
#' provided.
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
#'
discover_patterns <- function(data, cols = tidyselect::everything(),
                              group, type = "ngram", pattern,
                              len = 2:5, gap = 1:3, min_support = 0.01,
                              min_freq = 2, start, end, contains) {
  check_missing(data)
  data <- extract_data(data)
  group <- extract_group(data, group)
  cols <- get_cols(rlang::enquo(cols), data)
  last <- FALSE
  if (is.character(group) && length(group) == 1L) {
    if (group[1L] == "last_obs") {
      last <- TRUE
    } else {
      cols <- setdiff(cols, group)
    }
  }
  data <- prepare_sequence_data(data, cols)
  if (last) {
    data <- extract_last(data$sequences, data$alphabet)
    group <- data$group
  }
  discover_patterns_(
    sequences = data$sequences,
    alphabet = data$alphabet,
    group = group,
    type = type,
    pattern = pattern,
    len = len,
    gap = gap,
    min_support = min_support,
    min_freq = min_freq,
    start = start,
    end = end,
    contains = contains
  )
}

discover_patterns_ <- function(sequences, alphabet, group, type, pattern, len,
                               gap, min_support, min_freq, start, end,
                               contains) {
  m <- ncol(sequences)
  type <- check_match(type, c("ngram", "gapped", "repeated"))
  check_range(len, scalar = FALSE, type = "integer", min = 1L, max = m)
  check_range(gap, scalar = FALSE, type = "integer", min = 1L, max = m - 2L)
  check_range(min_support)
  check_values(min_freq)
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
  patterns <- format_patterns(patterns)
  patterns |>
    process_patterns(group) |>
    filter_patterns(min_support, min_freq, start, end, contains) |>
    pattern_proportions() |>
    pattern_lift(support = support) |>
    structure(
      patterns = patterns,
      group = group,
      class = c("patterns", "tbl_df", "tbl", "data.frame")
    )
}

process_patterns <- function(x, group) {
  out <- vector(mode = "list", length = length(x))
  groups <- NULL
  has_group <- FALSE
  if (!is.null(group)) {
    groups <- unique(group)
    num_grp <- as.integer(factor(group)) - 1L
    idx_grp <- lapply(groups, function(y) which(group == y))
    g <- length(groups)
    has_group <- TRUE
  }
  for (i in seq_along(out)) {
    pat_mat <- x[[i]]$patterns
    pat_len <- x[[i]]$length
    pat_u <- x[[i]]$unique
    u <- length(pat_u)
    if (u == 0) {
      next
    }
    n <- nrow(pat_mat)
    p <- ncol(pat_mat)
    pat_count <- .colSums(pat_mat > 0, m = n, n = p)
    pat_freq <- .colSums(pat_mat, m = n, n = p)
    if (has_group) {
      pat_count_grp <- matrix(0L, u, g)
      for (j in seq_len(g)) {
        idx <- idx_grp[[j]]
        m <- length(idx)
        pat_count_grp[, j] <- .colSums(pat_mat[idx, ] > 0, m = m, n = p)
      }
      pat_sep <- NULL
      if (g == 2L) {
        pat_sep <- .colSums(
          1L * (pat_mat > 0) == num_grp,
          m = n,
          n = p
        )
      }
    }
    tmp <- data.frame(
      pattern = pat_u,
      length = pat_len,
      frequency = pat_freq,
      count = pat_count
    )
    if (has_group) {
      colnames(pat_count_grp) <- paste0("count_", groups)
      tmp <- cbind(tmp, as.data.frame(pat_count_grp))
      chisq <- chisq_test(x = pat_count_grp, count = pat_count)
      tmp$chisq <- chisq$statistic
      tmp$p_value <- chisq$p_value
      if (g == 2L) {
        tmp$separation <- 2 * abs(0.5 - pat_sep / n)
      }
    }
    tmp$support <- tmp$count / n
    out[[i]] <- tmp
  }
  proto <- data.frame(
    pattern = character(0L),
    length = integer(0L),
    frequency = integer(0L),
    count = integer(0L),
    support = numeric(0L)
  )
  if (has_group) {
    proto_grp <- matrix(NA, nrow = 0L, ncol = g)
    colnames(proto_grp) <- paste0("count_", groups)
    proto <- cbind(proto, as.data.frame(proto_grp))
  }
  out <- dplyr::bind_rows(out, proto) |>
    dplyr::arrange(dplyr::desc(!!rlang::sym("frequency"))) |>
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
          base::paste,
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
          base::paste,
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
          base::paste,
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
        base::paste,
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

filter_patterns <- function(patterns, min_support, min_freq,
                            start, end, contains) {
  out <- patterns |>
    dplyr::filter(!!rlang::sym("support") >= min_support) |>
    dplyr::filter(!!rlang::sym("frequency") >= min_freq)
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

format_patterns <- function(x) {
  for (i in seq_along(x)) {
    pat_mat <- x[[i]]$patterns
    n <- nrow(pat_mat)
    p <- ncol(pat_mat)
    pat_vec <- c(pat_mat[nzchar(pat_mat)])
    if (length(pat_vec) == 0) {
      x[[i]]$patterns <- matrix(0L, nrow = n, ncol = 0L)
      x[[i]]$unique <- character(0L)
      next
    }
    pat_u <- unique(pat_vec)
    pat_out <- matrix(0L, n, length(pat_u))
    for (j in seq_len(p)) {
      pats <- pat_mat[, j]
      valid <- which(nzchar(pats))
      if (length(valid) > 0) {
        pats <- pats[valid]
        idx <- cbind(valid, match(pats, pat_u))
        pat_out[idx] <- pat_out[idx] + 1L
      }
    }
    x[[i]]$patterns <- pat_out
    x[[i]]$unique <- pat_u
  }
  x
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
      proportion = !!rlang::sym("frequency") / sum(!!rlang::sym("frequency"))
    ) |>
    dplyr::relocate(
      !!rlang::sym("proportion"),
      .after = !!rlang::sym("frequency")
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
  patterns |>
    dplyr::relocate(
      !!rlang::sym("lift"),
      .after = !!rlang::sym("support")
    )
}

# Modeling ----------------------------------------------------------------


#' Analyze Pattern-Outcome Relationships
#'
#' Fit a (mixed) logistic regression model to the data using sequence patterns
#' as predictors. The patterns are first determined using [discover_patterns()]
#' and used as predictors as specified by the user (either count or presence).
#' The user can further select the maximum number of patterns to use and how
#' to prioritize the patterns (separation, frequency, proportion of sequences,
#' etc.). Additional covariates and random effects can also be included in the
#' model. The logistic model is fitted using [stats::glm()] or by
#' [lme4::glmer()] in the case of a mixed model.
#'
#' @export
#' @inheritParams convert
#' @param group \[`expression`]\cr An optional tidy selection of columns that
#'   define the grouping factors. If provided, group-specific parameters
#'   can be specified via `re_formula`. The default value `NULL` disables
#'   grouping and a fixed-effects model is fitted instead.
#' @param outcome \[`character(1)`]\cr Outcome variable specification. The
#'   option `"last_obs"` (default) assumes that the last non-missing observation
#'   of each sequence specifies the outcome. This can also be a column name of
#'   `data` that contains the outcome variable information.
#' @param n \[`integer(1)`]\cr Maximum number of patterns to include in the
#'   model as covariates. The default is `10`.
#' @param type \[`character(1)`]\cr Specifies how the patterns are included
#'   as predictors. The option `"count"` uses the number of times a specific
#'   pattern occurs in a sequence, whereas option `"presence"` (the default)
#'   defines a binary predictors that attains the value `1` if the pattern is
#'   present in the sequence and `0` otherwise.
#' @param priority \[`character(1)`]\cr Priority of pattern inclusion in the
#'   model. The available options are:
#'
#'   * `"separation"`: Best separation of the `outcome` classes. For
#'     `type = "count"`, the chi-square statistic is used. For
#'     `type = "presence"`, the proportion of 1-1 and 0-0 matches between
#'     the pattern and the `outcome` is used.
#'   * `"frequency"`: Most frequent patterns.
#'   * `"support"`: Most common patterns.
#'
#' @param formula \[`formula`]\cr Formula specification for the fixed effects
#'   of the non-pattern covariates. The default is `~ 1`, adding no covariates.
#' @param re_formula \[`formula`]\cr Formula specification of the random
#'   effects when `group` is provided. By default, a random intercept is added
#'   for each grouping variable in `group`.
#' @param ... Additional arguments passed to [discover_patterns()] and the
#'   model-fitting function ([stats::glm()] or [lme4::glmer()]).
#' @return Either a `glm` or a `glmerMod` object depending on whether
#'   random effects were included.
#' @examples
#' # TODO
#' 1 + 1
#'
analyze_outcome <- function(data, cols = tidyselect::everything(),
                            group = NULL, outcome = "last_obs",
                            n = 10, type = "presence",
                            priority = "separation", formula = ~1,
                            re_formula, ...) {
  check_missing(data)
  check_values(n, strict = TRUE)
  type <- check_match(type, c("presence", "count"))
  priority <- check_match(priority, c("separation", "frequency", "support"))
  data <- extract_data(data)
  group <- get_cols(rlang::enquo(group), data)
  cols <- get_cols(rlang::enquo(cols), data) |>
    setdiff(c(group, outcome))
  mixed <- !is.null(group)
  seqdata <- prepare_sequence_data(data, cols)
  if (identical(outcome, "last_obs")) {
    seqdata <- extract_last(seqdata$sequences, seqdata$alphabet)
    outcome <- seqdata$group
    data$.outcome <- factor(outcome)
    response <- ".outcome"
  } else {
    outcome <- as.character(outcome)
    stopifnot_(
      length(outcome) == 1L && outcome %in% names(data),
      "Argument {.arg outcome} must be a column name of {.arg data}."
    )
    response <- outcome
    data[[outcome]] <- factor(data[[outcome]])
    outcome <- data[[outcome]]
  }
  stopifnot_(
    n_unique(outcome) == 2L,
    "Argument {.arg outcome} must specify an outcome variable with two classes."
  )
  dots <- list(...)
  dp_args_names <- c(
    "type", "pattern", "len",
    "gap", "min_support", "min_freq",
    "start", "end", "contains"
  )
  dots_names <- names(dots)
  fit_args <- dots[setdiff(dots_names, dp_args_names)]
  dp_args <- dots[intersect(dots_names, dp_args_names)]
  dp_args$sequences <- seqdata$sequences
  dp_args$alphabet <- seqdata$alphabet
  dp_args$group <- outcome
  dp_args$type <- dp_args$type %||% "ngram"
  dp_args$len <- dp_args$len %||% 1:2
  dp_args$gap <- dp_args$gap %||% 1
  dp_args$min_support <- dp_args$min_support %||% 0.05
  dp_args$min_freq <- dp_args$min_freq %||% 10L
  disc <- do.call(
    getFromNamespace("discover_patterns_", "codyna"),
    dp_args
  ) |>
    dplyr::arrange(dplyr::desc(!!rlang::sym(priority))) |>
    dplyr::slice_head(n = n)
  filtered <- disc$pattern
  patterns <- attr(disc, "patterns") |>
    lapply(
      function(x) {
        out <- x$patterns
        colnames(out) <- x$unique
        out[, x$unique %in% filtered]
      }
    ) |>
    do.call(what = base::cbind) |>
    as.data.frame()
  if (type == "presence") {
    patterns <- patterns |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),
          function(x) 0L + 1L * (x > 0)
        )
      )
  }
  pat_formula <- paste0("`", filtered, "`") |>
    paste0(collapse = " + ") |>
    paste("~ ", ... = _) |>
    stats::as.formula()
  intercept <- attr(terms(formula), "intercept") == 1
  terms <- union(
    attr(stats::terms(formula), "term.labels"),
    attr(stats::terms(pat_formula), "term.labels")
  )
  formula <- stats::reformulate(
    termlabels = terms,
    response = response,
    intercept = intercept
  )
  fit_fun <- stats::glm
  if (mixed) {
    re_formula <- re_formula %m%
      paste0("(1 | ", group, ")", collapse = " + ") |>
      paste("~ . + ", ... = _) |>
      stats::as.formula()
    formula <- stats::update(formula, re_formula)
    fit_fun <- lme4::glmer
  }
  fit_args$formula <- formula
  fit_args$family <- stats::binomial()
  fit_args$data <- cbind(data, patterns)
  fit <- try_(do.call(fit_fun, fit_args))
  stopifnot_(
    !inherits(fit, "try-error"),
    "Model fitting failed"
  )
  if (mixed) {
    fit@call$data <- NULL
  } else {
    fit$call <- NULL
  }
  fit
}
