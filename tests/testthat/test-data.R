test_that("sequence data can be converted to other formats", {
  convert(mock_sequence, format = "frequency") |>
    expect_error(NA)
  convert(mock_sequence, format = "onehot") |>
    expect_error(NA)
  convert(mock_sequence, format = "edgelist") |>
    expect_error(NA)
  convert(mock_sequence, format = "reverse") |>
    expect_error(NA)
})

test_that("tna objects can be converted to sequence data", {
  cols <- colnames(mock_sequence_num)
  data <- extract_data(mock_tna) |>
    prepare_sequence_data(cols) |>
    expect_error(NA)
  expect_equal(data$sequences, mock_sequence_num, ignore_attr = TRUE)
  expect_equal(data$alphabet, c("A", "B", "C"))
})

test_that("data extraction fails from a matrix without colnames", {
  mat <- as.matrix(mock_sequence)
  dimnames(mat) <- NULL
  extract_data(mat) |>
    expect_error(
      "Argument `data` must have column names when a <matrix> is provided"
    )
})
