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
  patterns <- discover_patterns(mock_sequence, contains = "C|B")
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
