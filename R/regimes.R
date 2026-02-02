#' Regime Detection for Time Series Data
#'
#' Detects regime changes in time series data using multiple methods including
#' cumulative peaks, changepoint detection, variance shifts,
#' threshold analysis, gradient changes, and entropy analysis.
#'
#' @export
#' @param data \[`ts`, `numeric()`]\cr Univariate time series data.
#' @param method \[`character(1)`]\cr Detection method.
#'   The available options are:
#'
#'   * `"cumulative_peaks"`: Detects cumulative complexity peaks using Z-tests.
#'   * `"changepoint"`: Change point detection (multi-window mean-shift test).
#'   * `"threshold"`: Adaptive quartile-based regime classification.
#'   * `"variance_shift"`: Detects changes in variance patterns.
#'   * `"slope"`: Detects changes in local slope (rolling linear models).
#'   * `"entropy"`: Detects changes in the Shannon entropy of the
#'     complexity series, calculated in rolling windows.
#'   * `"smart"` (default): Combines gradient, peaks, and changepoint methods.
#'   * `"all"`: Applies all individual methods listed above and uses
#'     ensemble voting.
#'
#' @param sensitivity \[`character(1)`]\cr Detection sensitivity level.
#'   The available options are: `"low"`, `"medium"`, `"high"`.
#'   The default is `"medium"`. This affects thresholds and window sizes within
#'   the detection methods.
#' @param min_change \[`integer(1)`]\cr Minimum number of observations between
#'   changes. If missing (default), the value is determined automatically
#'   (typically 10% of observations, minimum of 10).
#' @param window \[`integer(1)`]\cr base window size for rolling calculations.
#'   This is further adjusted by `sensitivity`. The default is `10`.
#' @param peak \[`numeric(1)`]\cr Base z-score threshold for individual peak
#'   detection with the `"cumulative_peaks"` method.
#'   Adjusted by `sensitivity`. The default is `2.0`.
#' @param cumulative \[`numeric(1)`]\cr A value between 0 and 1 that defines
#'   the base proportion threshold for identifying cumulative peak regions.
#'   Adjusted by `sensitivity`. The default is `0.6`.
#'
#' @return An object of class `regimes` which is a `tibble` containing
#'   the following columns:
#'
#'   `value`: Original time series data.
#'   `time`: Original time points.
#'   `change`: A `logical` vector indicating regime changes.
#'   `id`: An `integer` regime identifier.
#'   `type`: Type of change detected by the method.
#'   `magnitude`: Magnitude of the change (method-specific interpretation)
#'   `confidence`: Confidence in the detection
#'     (method-specific interpretation, typically between `0` and `1`, or `NA`)
#'   `stability`: Categorical stability: `"Stable"`, `"Transitional"`, and
#'     `"Unstable"`.
#'   `score`: A `numeric` stability score between `0` and `1`.
#'
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' regimes <- detect_regimes(
#'   data = ts_data,
#'   method = "threshold",
#'   sensitivity = "medium"
#' )
#'
detect_regimes <- function(data, method = "smart", sensitivity = "medium",
                           min_change, window = 10, peak = 2.0,
                           cumulative = 0.6) {
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  method <- check_match(method, detection_methods)
  sensitivity <- check_match(sensitivity, c("low", "medium", "high"))
  n <- length(values)
  check_range(window, type = "integer", min = 2L, max = n - 1)
  check_range(cumulative)
  check_values(peak, type = "numeric")
  min_change <- ifelse_(missing(min_change), max(10, floor(n * 0.10)))
  min_change <- max(1, floor(min_change))
  params <- default_detection_parameters(
    sensitivity,
    window,
    peak,
    cumulative,
    min_change
  )
  result <- do.call(
    what = getFromNamespace(paste0("detect_", method), "codyna"),
    args = list(values = values, time = time, params = params)
  )
  orig <- data.frame(
    value = values,
    time = time
  )
  regimes <- cbind(orig, result)
  stb <- regime_stability(regimes, min_change)
  regimes$stability <- stb$stability
  regimes$score <- stb$score
  structure(
    tibble::as_tibble(regimes),
    class = c("regimes", "tbl_df", "tbl", "data.frame")
  )
}

default_detection_parameters <- function(sensitivity, window, peak,
                                         cumulative, min_change) {
  params <- list(
    window = max(3, floor(window)),
    peak = peak,
    cumulative = cumulative,
    effect = 0.5,
    variance_ratio = 1.5,
    slope_signif = 0.2,
    slope_factor = 0.005,
    entropy_bins = 10,
    entropy_rel_change = 0.2,
    sensitivity = sensitivity,
    min_change = min_change
  )
  if (sensitivity == "high") {
    params$window <- max(3, floor(params$window * 0.7))
    params$peak <- params$peak * 0.6
    params$cumulative <- params$cumulative * 0.5
    params$effect <- params$effect * 0.6
    params$variance_ratio <- 1 + (params$variance_ratio - 1) * 0.6
    params$slope_signif <- params$slope_signif * 0.7
    params$slope_factor <- params$slope_factor * 0.7
    params$entropy_bins <- max(5, floor(params$entropy_bins * 0.7))
    params$entropy_rel_change <- params$entropy_rel_change * 0.6
  } else if (sensitivity == "low") {
    params$window <- max(5, floor(params$window_size * 1.3))
    params$peak <- params$peak * 1.4
    params$cumulative <- min(1, params$cumulative * 1.3)
    params$effect <- params$effect * 1.4
    params$variance_ratio <- 1 + (params$variance - 1) * 1.4
    params$slope_signif <- params$slope_signif * 1.3
    params$slope_factor <- params$slope_factor * 1.3
    params$entropy_bins <- floor(params$entropy_bins * 1.3)
    params$entropy_rel_change <- params$entropy_rel_change * 1.4
  }
  params$window <- ifelse_(
    params$window %% 2L == 0L,
    params$window + 1L,
    params$window
  )
  params$entropy_nbins <- max(2, params$entropy_bins)
  params
}

min_change_constraint <- function(change, min_change) {
  if (min_change <= 1 || !any(change)) {
    return(change)
  }
  change_idx <- which(change)
  n <- length(change_idx)
  new_change <- rep(FALSE, length(change))
  last_idx <- change_idx[1L]
  if (length(change_idx) > 1) {
    for (i in seq(2L, n)) {
      current_idx <- change_idx[i]
      if ((current_idx - last_idx) >= min_change) {
        new_change[current_idx] <- TRUE
        last_idx <- current_idx
      }
    }
  }
  new_change
}

regime_stability <- function(regimes, min_change) {
  n <- nrow(regimes)
  score <- rep(NA_real_, n)
  stability <- rep("Unstable", n)
  u_id <- setdiff(unique(regimes$id), 1)
  for (id in u_id) {
    idx <- which(regimes$id == id)
    len <- length(idx)
    for (j in seq_len(len)) {
      k <- idx[j]
      denom <- max(1, floor(len / 2.0))
      dist <- min(j - 1, len - j) / denom
      boundary <- min(1, dist)
      duration <- min(1, len / max(min_change, 10))
      score[k] <- (boundary + duration) / 2.0
    }
  }
  stability[score >= 0.35] <- "Transitional"
  stability[score >= 0.75] <- "Stable"
  stability[regimes$id == 1] <- "Initial"
  list(stability = stability, score = score)
}

# Detection methods -------------------------------------------------------

detection_methods <- c(
  "all",
  "cumulative_peaks",
  "changepoints",
  "entropy",
  "slope",
  "smart",
  "threshold",
  "variance_shift"
)

detect_all <- function(values, time, params) {
  methods <- setdiff(detection_methods, c("all", "smart"))
  results <- list()
  n <- length(values)
  args <- list(values = values, time = time, params = params)
  for (method in methods) {
    results[[method]] <- do.call(
      what = getFromNamespace(paste0("detect_", method), "codyna"),
      args = args
    )
  }
  changes <- lapply(results, "[[", "change")
  vote_matrix <- do.call(base::cbind, changes)
  ensemble_vote_threshold_prop <- 0.3
  vote_prop <- (params$sensitivity == "medium") * 0.3 +
    (params$sensitivity == "high") * 0.15 +
    (params$sensitivity == "low") * 0.5
  n_voted <- ncol(vote_matrix)
  vote_sums <- rowSums(vote_matrix, na.rm = TRUE)
  min_votes <- max(1, floor(ensemble_vote_threshold_prop * n_voted))
  change <- vote_sums >= min_votes
  change <- min_change_constraint(change, params$min_change)
  type <- rep("none", n)
  magnitude <- vote_sums / max(1, n_voted)
  confidence <- vote_sums / max(1, n_voted)
  type[change] <- paste0(
    "ensemble_vote_",
    round(magnitude[change] * 100),
    "%"
  )
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_cumulative_peaks <- function(values, time, params) {
  n <- length(values)
  window <- params$window
  peak <- params$peak
  cumulative <- params$cumulative
  peak_indicators <- rep(FALSE, n)
  z_scores <- rep(NA, n)
  half_window <- window %/% 2L
  for (i in seq_len(n)) {
    val <- values[i]
    if (is.na(val)) {
      next
    }
    start <- max(1L, i - half_window)
    end <- min(n, i + half_window)
    idx <- c(start:(i - 1L), (i + 1L):end)
    window_data <- values[idx]
    window_data <- window[!is.na(window_data)]
    if (length(window_data) >= 3) {
      window_mean <- mean(window_data, na.rm = TRUE)
      window_sd <- stats::sd(window_data, na.rm = TRUE)
      if (!is.na(window_sd)) {
        z <- (val - window_mean) / window_sd
        z_scores[i] <- z
        if (!is.na(z) && is.finite(z) && z > peak) {
          peak_indicators[i] <- TRUE
        }
      }
    }
  }
  peak_region <- rep(FALSE, n)
  peak_proportion <- rollmean(peak_indicators, window, align = "center")
  peak_valid <- !is.na(peak_proportion)
  peak_region[peak_valid] <- peak_proportion[peak_valid] >= cumulative
  change <- rep(FALSE, n)
  change[seq(2L, n)] <- diff(peak_region) != 0
  change <- min_change_constraint(change, params$min_change)
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- rep(0, n)
  change_idx <- which(change)
  for (i in change_idx) {
    if (i == 1) {
      next
    }
    type[i] <- ifelse(
      peak_region[i],
      "peak_region_onset",
      "peak_region_offset"
    )
    magnitude[i] <- ifelse_(
      is.na(z_scores[i]),
      0,
      abs(z_scores[i])
    )
    confidence[i] <- min(
      1,
      ifelse_(
        !is.na(z_scores[i]) && peak > 0,
        abs(z_scores[i]) / peak,
        0.5
      )
    )
  }
  magnitude[change & is.na(magnitude)] <- 0
  confidence[change & (is.na(confidence) | confidence == 0)] <- 0.5
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_changepoints <- function(values, time, params) {
  n <- length(values)
  effect <- params$effect
  window <- params$window
  change_points <- integer(0)
  window_mult <- c(0.75, 1.0, 1.5)
  for (mult in window_mult) {
    segment_len <- max(3, floor(window * mult))
    if (n < segment_len * 2 + 1) {
      next
    }
    for (i in (segment_len + 1):(n - segment_len)) {
      segment_before <- values[(i - segment_len):(i - 1L)]
      segment_after <- values[i:(i + segment_len - 1L)]
      segment_before <- segment_before[!is.na(segment_before)]
      segment_after <- segment_after[!is.na(segment_after)]
      if (length(segment_before) >= 2 && length(segment_after) >= 2) {
        mean_b <- mean(segment_before, na.rm = TRUE)
        mean_a <- mean(segment_after, na.rm = TRUE)
        sd_b <- stats::sd(segment_before, na.rm = TRUE)
        sd_a <- stats::sd(segment_after, na.rm = TRUE)
        combined_sd <- stats::sd(
          c(segment_before, segment_after), na.rm = TRUE
        )
        effect_size <- abs(mean_a - mean_b) / combined_sd
        if (!is.na(effect_size) && effect_size > effect) {
          change_points <- c(change_points, i)
        }
      }
    }
  }
  change <- rep(FALSE, n)
  if (length(change_points) > 0) {
    unique_change_points <- sort(unique(change_points))
    tmp <- rep(FALSE, n)
    tmp[unique_change_points] <- TRUE
    change <- min_change_constraint(tmp, params$min_change)
  }
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- rep(0, n)
  change_idx <- which(change)
  m <- max(3, floor(window * 0.5))
  sd_ref <- stats::sd(values, na.rm = TRUE)
  for (i in change_idx){
    if (i > m && i <= (n - m + 1)){
      mean_l <- mean(values[max(1, i - m):(i - 1)], na.rm = TRUE)
      mean_r <- mean(values[i:min(n, i + m - 1)], na.rm = TRUE)
      if (!is.na(mean_l) && !is.na(mean_r)){
        type[i] <- ifelse_(
          mean_r > mean_l,
          "mean_level_increase",
          "mean_level_decrease"
        )
        mag <- abs(mean_r - mean_l) / sd_ref
        magnitude[i] <- mag
        confidence[i] <- min(1, mag / effect)
      } else {
        type[i] <- "mean_level_change"
        confidence[i] <- 0.5
      }
    } else {
      type[i] <- "mean_level_change_edge"
      confidence[i] <- 0.5
    }
  }
  confidence[change & (is.na(confidence) | confidence == 0)] <- 0.5
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_entropy <- function(values, time, params) {
  n <- length(values)
  window <- params$window
  bins <- params$entropy_bins
  rel <- params$entropy_rel_change
  ent <- roll(entropy, values, window = window, align = "center", bins = bins)
  change <- rep(FALSE, n)
  for (i in 2:n) {
    if (!is.na(ent[i]) && !is.na(ent[i - 1])) {
      rel_diff <- abs(ent[i] - ent[i - 1]) / ent[i - 1]
      change[i] <- rel_diff > rel
    }
  }
  change <- min_change_constraint(change, params$min_change)
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- rep(0, n)
  change_idx <- which(change)
  for (i in change_idx){
    if (i > 1 && !is.na(ent[i]) && !is.na(ent[i - 1])) {
      type[i] <- ifelse_(
        ent[i] > ent[i - 1],
        "entropy_increase",
        "entropy_decrease"
      )
      abs_diff_ent <- abs(ent[i] - ent[i - 1])
      magnitude[i] <- abs_diff_ent
      rel_diff <- abs_diff_ent / ent[i - 1]
      confidence[i] <- min(1, rel_diff / rel)
    } else {
      type[i] <- "entropy_change"
      confidence[i] <- 0.5
    }
  }
  confidence[change & (is.na(confidence) | confidence == 0)] <- 0.5
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_slope <- function(values, time, params) {
  n <- length(values)
  window <- params$window
  slope_signif <- params$slope_signif
  slope_factor <- params$slope_factor
  slope <- rep(NA, n)
  r_squared <- rep(NA, n)
  half_window <- window %/% 2
  for (i in 1:n) {
    start <- max(1, i - half_window)
    end <- min(n, i + half_window)
    w <- start:end
    if (length(w) < 3) {
      slope[i] <- 0
      r_squared[i] <- 0
      next
    }
    y <- values[w]
    x <- time[w]
    obs <- !is.na(y)
    y <- y[obs]
    x <- x[obs]
    if (length(y) >= 3 && length(unique(y)) > 1) {
      fit <- stats::lm(y ~ x)
      if (!is.null(fit) && length(stats::coef(fit)) == 2
          && !any(is.na(stats::coef(fit)))) {
        slope[i] <- stats::coef(fit)[2L]
        r_squared[i] <- summary(fit)$r.squared
      } else {
        slope[i] <- 0
        r_squared[i] <- 0
      }
    } else {
      slope[i] <- 0
      r_squared[i] <- 0
    }
  }
  change <- rep(FALSE, n)
  sd_ref <- stats::sd(values, na.rm = TRUE)
  slope_thr <- slope_factor * sd_ref
  for (i in 2L:n) {
    changed <- FALSE
    if (is.na(slope[i]) || is.na(slope[i - 1]) ||
        is.na(r_squared[i]) || is.na(r_squared[i - 1])) {
      next
    }
    if (r_squared[i] > slope_signif && r_squared[i - 1] > slope_signif &&
        sign(slope[i]) != sign(slope[i-1]) &&
        abs(slope[i]) > slope_thr && abs(slope[i - 1]) > slope_thr) {
      changed <- TRUE
    }
    if (abs(slope[i] - slope[i - 1]) > (2 * slope_thr)) {
      changed <- TRUE
    }
    sig <- r_squared[i] > slope_signif && abs(slope[i]) > slope_thr
    sig_prev <- r_squared[i - 1] > slope_signif && abs(slope[i - 1]) > slope_thr
    if (sig != sig_prev) {
      changed <- TRUE
    }
    change[i] <- changed
  }
  change <- min_change_constraint(change, params$min_change)
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- rep(0, n)
  change_idx <- which(change)
  for (i in change_idx) {
    s <- slope[i]
    s_prev <- slope[max(1, i - 1)]
    rsq <- r_squared[i]
    if (!is.na(s) && !is.na(s_prev) && !is.na(rsq)){
      if (rsq > slope_signif * 0.75 && abs(s) > slope_thr * 0.5) {
        if (abs(s) < slope_thr * 0.75) type[i] <- "trend_flattening"
        else type[i] <- ifelse_(s > 0, "trend_increasing", "trend_decreasing")
      } else {
        type[i] <- "trend_shift_subtle"
      }
      magnitude[i] <- abs(s - s_prev) / sd_ref
      conf_val <- rsq / slope_signif + abs(s - s_prev) / slope_thr
      confidence[i] <- min(1, max(0, conf_val, na.rm = TRUE), na.rm = TRUE)
    } else {
      type[i] <- "trend_change"
      confidence[i] <- 0.3
    }
  }
  confidence[change & (is.na(confidence) | confidence == 0)] <- 0.3
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_smart <- function(values, time, params) {
  n <- length(values)
  res_slope <- detect_slope(values, time, params)
  res_peaks <- detect_cumulative_peaks(values, time, params)
  res_change <- detect_changepoints(values, time, params)
  w_slope <- 0.4
  w_peaks <- 0.3
  w_change <- 0.3
  combined_score <- res_slope$change * w_slope +
    res_peaks$change * w_peaks + res_change$change * w_change
  consensus_prop <- (params$sensitivity == "medium") * 0.35 +
    (params$sensitivity == "high") * 0.25 +
    (params$sensitivity == "low") * 0.5
  max_score <- w_slope + w_peaks + w_change
  change <- combined_score >= (consensus_prop * max_score)
  change <- min_change_constraint(change, params$min_change)
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- pmin(1, pmax(0, combined_score / max_score))
  change_idx <- which(change)
  for (i in change_idx) {
    type_set <- FALSE
    if (!is.null(res_slope) && res_slope$change[i]
        && res_slope$type[i] != "none") {
      type[i] <- res_slope$type[i]
      magnitude[i] <- res_slope$magnitude[i]
      type_set <- TRUE
    } else if (!is.null(res_peaks) && res_peaks$change[i]
               && res_peaks$type[i] != "none") {
      type[i] <- res_peaks$type[i]
      magnitude[i] <- res_peaks$magnitude[i]
      type_set <- TRUE
    } else if (!is.null(res_change) && res_change$change[i]
               && res_change$type[i] != "none") {
      type[i] <- res_change$type[i]
      magnitude[i] <- res_change$magnitude[i]
      type_set <- TRUE
    }
    if (!type_set) {
      type[i] <- "smart_combo_general"
      if (magnitude[i] == 0 && !is.na(confidence[i])) {
        magnitude[i] <- confidence[i]
      }
    }
  }
  confidence[change & (is.na(confidence) | confidence == 0)] <- 0.5
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_threshold  <- function(values, time, params) {
  n <- length(values)
  probs <- c(0, 1/3, 2/3, 1.0)
  q <- unique(stats::quantile(values, probs = probs, na.rm = TRUE))
  nq <- length(q)
  if (nq < 2L) {
    warning_("Could not define distinct quantiles for threshold method.")
    regimes_raw <- rep(1L, n)
  } else {
    regimes_raw <- as.integer(
      cut(values, breaks = q, include.lowest = TRUE, labels = FALSE)
    )
  }
  regimes_raw[is.na(regimes_raw)] <- 1L
  smoothed_regimes <- regimes_raw
  m <- max(3, params$window %/% 2)
  if (m %% 2 == 0) {
    m <- m + 1
  }
  if (n >= m) {
    mode_roll <- roll(stat_mode, regimes_raw, window = m, align = "center")
    valid <- !is.na(mode_roll)
    smoothed_regimes[valid] <- mode_roll[valid]
  }
  smoothed_regimes[is.na(smoothed_regimes)] <- 1L
  change <- rep(FALSE, n)
  change[2:n] <- diff(smoothed_regimes) != 0
  changes <- min_change_constraint(change, params$min_change)
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- rep(NA, n)
  num_levels <- nq - 1
  if (is.na(num_levels) || num_levels <= 0) {
    num_levels <- 1
  }
  lab <- paste0("Level ", seq_len(num_levels))
  if (num_levels == 3) {
    lab <- c("low", "medium", "high")
  }
  change_idx <- which(change)
  for (i in change_idx) {
    if (i == 1) {
      next
    }
    from <- smoothed_regimes[i - 1]
    to <- smoothed_regimes[i]
    if(!is.na(from) && !is.na(to) && from >= 1 && from <= num_levels &&
       to >= 1 && to <= num_levels) {
      type[i] <- paste0("threshold_", lab[from], "_to_", lab[to])
      magnitude[i] <- abs(to - from) / max(1, nq)
    } else {
      type[i] <- "threshold_level_change"
      magnitude[i] <- 0.5
    }
  }
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}

detect_variance_shift <- function(values, time, params) {
  n <- length(values)
  window <- params$window
  var_ratio <- params$variance_ratio
  if (n < window * 2 + 1) {
    warning_(
      "Not enough data for {.var variance_shift} method
       with window size {window}."
    )
    return(
      list(
        change = rep(FALSE, n),
        id = rep(1L, n),
        type = rep("var_insufficient_data", n),
        magnitude = rep(0, n),
        confidence = rep(0, n)
      )
    )
  }
  var_ratio_log <- rep(NA, n)
  for (i in (window + 1):(n - window)) {
    segment_before <- values[(i - window):(i - 1)]
    segment_after <- values[i:(i + window - 1)]
    segment_before <- segment_before[!is.na(segment_before)]
    segment_after <- segment_after[!is.na(segment_after)]
    if (length(segment_before) >= 2 && length(segment_after) >= 2) {
      var_before <- stats::var(segment_before, na.rm = TRUE)
      var_after <- stats::var(segment_after, na.rm = TRUE)
      if (!is.na(var_before) && !is.na(var_after)) {
        var_ratio_log[i] <- log(var_after / var_before)
      }
    }
  }
  change <- rep(FALSE, n)
  potential <- which(
    !is.na(var_ratio_log) & abs(var_ratio_log) > log(var_ratio)
  )
  if(length(potential) > 0){
    tmp <- rep(FALSE, n)
    tmp[potential] <- TRUE
    change <- min_change_constraint(tmp, params$min_change)
  }
  type <- rep("none", n)
  magnitude <- rep(0, n)
  confidence <- rep(0, n)
  change_idx <- which(change)
  for (i in change_idx) {
    r <- var_ratio_log[i]
    if (!is.na(r)){
      type[i] <- ifelse_(r > 0, "variance_increase", "variance_decrease")
      magnitude[i] <- abs(r)
      confidence[i] <- min(1, abs(r) / log(var_ratio))
    } else {
      type[i] <- "variance_change"
      confidence[i] <- 0.5
    }
  }
  confidence[change & (is.na(confidence) | confidence == 0)] <- 0.5
  list(
    change = change,
    id = 1L + cumsum(change),
    type = type,
    magnitude = magnitude,
    confidence = confidence
  )
}
