test_that("sequence indices can be computed", {
  sequence_indices(mock_sequence) |>
    expect_error(NA)
})

test_that("sequence indices supports different input formats", {
  sequence_indices(engagement[1:100, ]) |>
    expect_error(NA)
  sequence_indices(mock_tna) |>
    expect_error(NA)
  sequence_indices(mock_group_tna) |>
    expect_error(NA)
  sequence_indices(mock_tna_data) |>
    expect_error(NA)
})

test_that("favorable states can be specified", {
  sequence_indices(mock_sequence, favorable = "A") |>
    expect_error(NA)
})
