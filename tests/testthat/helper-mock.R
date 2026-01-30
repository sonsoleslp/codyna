set.seed(0)

mock_sequence <- data.frame(
  T1 = c("A", "B", "C", "A", "A"),
  T2 = c("A", "C", "B", "B", "B"),
  T3 = c("B", "C", "A", "C", "A"),
  T4 = c("B", "A", "C", "A", "B"),
  T5 = c("C", "A", "B", "B", "A"),
  T6 = c("C", "C", "A", "A", "C")
)

mock_ts <- stats::arima.sim(
  list(order = c(2, 1, 0), ar = c(0.5, 0.2)),
  n = 100
)

mock_sequence_num <- c(A = 1, B = 2, C = 3)[as.matrix(mock_sequence)]
dim(mock_sequence_num) <- dim(mock_sequence)
colnames(mock_sequence_num) <- paste0("T", 1:6)

mock_tna <- structure(
  list(
    weights = 0,
    data = structure(
      mock_sequence_num,
      alphabet = c("A", "B", "C")
    )
  ),
  class = "tna"
)
