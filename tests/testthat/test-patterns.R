test_that("n-gram patterns can be discovered", {
  discover_patterns(mock_sequence, type = "ngram") |>
    expect_error(NA)
})

test_that("gapped pattern can be discovered", {
  discover_patterns(mock_sequence, type = "gapped") |>
    expect_error(NA)
})

test_that("repeated pattern can be discovered", {
  discover_patterns(mock_sequence, type = "repeated") |>
    expect_error(NA)
})

test_that("custom patterns can be discovered", {
  discover_patterns(mock_sequence, pattern = "A->*") |>
    expect_error(NA)
  discover_patterns(mock_sequence, pattern = "*->B") |>
    expect_error(NA)
})

test_that("patterns can be filtered based on states", {
  patterns <- discover_patterns(mock_sequence, start = "B", end = "C")
  expect_true("B->A->C" %in% patterns$pattern)
  expect_true("B->C->C" %in% patterns$pattern)
  patterns <- discover_patterns(mock_sequence, contain = "C|B")
  expect_true("A->C" %in% patterns$pattern)
  expect_true("A->B" %in% patterns$pattern)
  expect_true("C->C" %in% patterns$pattern)
  expect_false("A->A" %in% patterns$pattern)
})

test_that("patterns can be filtered based on frequency and support", {
  patterns <- discover_patterns(mock_sequence, min_freq = 3)
  expect_true(all(patterns$frequency >= 3))
  patterns <- discover_patterns(mock_sequence, min_support = 0.2)
  expect_true(all(patterns$support >= 0.2))
})

test_that("output has zero rows if no patterns are found", {
  patterns <- discover_patterns(mock_sequence, pattern = "A->******->B")
  expect_true(nrow(patterns) == 0L)
  patterns <- discover_patterns(mock_sequence, min_freq = 6)
  expect_true(nrow(patterns) == 0L)
})

test_that("pattern discovery supports different input formats", {
  discover_patterns(engagement[1:100, ]) |>
    expect_error(NA)
  discover_patterns(mock_tna) |>
    expect_error(NA)
  discover_patterns(mock_group_tna) |>
    expect_error(NA)
  discover_patterns(mock_tna_data) |>
    expect_error(NA)
})

test_that("pattern discovery supports last observation as outcome", {
  mock_sequence_out <- mock_sequence
  mock_sequence_out$T7 <- rep(c("D", "F"), length.out = nrow(mock_sequence))
  patterns <- discover_patterns(mock_sequence_out, outcome = "last_obs") |>
    expect_error(NA)
  all(c("count_D", "count_F", "chisq") %in% names(patterns)) |>
    expect_true()
})

test_that("discovered patterns can be printed", {
  patterns <- discover_patterns(mock_sequence)
  print(patterns) |>
    expect_error(NA)
})

test_that("discovered patterns can be plotted", {
  patterns1 <- discover_patterns(mock_sequence)
  plot(patterns1) |>
    expect_error(NA)
  out <- rep(1:2, length.out = nrow(mock_sequence))
  patterns2 <- discover_patterns(mock_sequence, outcome = out)
  plot(patterns2) |>
    expect_error(NA)
  plot(patterns2, group = "1", global = FALSE) |>
    expect_error(NA)
  plot(patterns2, group = "1", global = TRUE) |>
    expect_error(NA)
})

test_that("outcomes can be analyzed", {
  fit <- analyze_outcome(
    engagement[idx, ],
    outcome = rep(1:2, length.out = length(idx)),
    len = 1
  ) |>
    expect_error(NA)
})

test_that("outcome analysis supports different input formats", {
  out <-  rep(1:2, length.out = nrow(mock_sequence))
  analyze_outcome(mock_tna, outcome = out, len = 2) |>
    expect_error(NA)
  analyze_outcome(mock_group_tna, outcome = out, len = 2, mixed = FALSE) |>
    expect_error(NA)
  analyze_outcome(mock_tna_data, outcome = out, len = 2) |>
    expect_error(NA)
})

test_that("outcome analysis supports mixed effects", {
  out <-  rep(1:2, length.out = nrow(mock_sequence))
  analyze_outcome(mock_tna_data, outcome = out, group = "group", len = 2) |>
    suppressMessages() |> # singular
    expect_error(NA)
  analyze_outcome(mock_group_tna, outcome = out, len = 2) |>
    suppressMessages() |> # singular
    expect_error(NA)
})
