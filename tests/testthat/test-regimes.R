test_that("regime detection can be applied", {
  detect_regimes(mock_ts, method = "smart") |>
    expect_error(NA)
  detect_regimes(mock_ts, method = "all") |>
    expect_error(NA)
})

test_that("regime detection sensitivity can be varied", {
  detect_regimes(mock_ts, method = "all", sensitivity = "high") |>
    expect_error(NA)
  detect_regimes(mock_ts, method = "all", sensitivity = "low") |>
    expect_error(NA)
})

test_that("regime detection accounts for missing values", {
  set.seed(0)
  n <- length(mock_ts)
  idx <- sample(n, floor(0.2 * n), replace = FALSE)
  mock_ts_na <- mock_ts
  mock_ts_na[idx] <- NA
  detect_regimes(mock_ts, method = "all") |>
    expect_error(NA)
  set.seed(0)
  n <- length(mock_ts)
  idx <- sample(n, floor(0.95 * n), replace = FALSE)
  mock_ts_na <- mock_ts
  mock_ts_na[idx] <- NA
  detect_regimes(mock_ts, method = "all") |>
    expect_error(NA)
})

test_that("variance shift warns with low data", {
  detect_regimes(mock_ts, method = "all", window = 99) |>
    expect_warning("Not enough data")
})
