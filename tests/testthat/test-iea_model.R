test_that("iea_model works", {
  A <-  matrix(c(1, 1, 0,
                 1, 2, 2,
                 0, 2, 0),
                nrow = 3, ncol = 3)


  res <- iea_model(adj = A , type = 'graph', model = 'IEAS', K = 0, apx = FALSE)
  expect_equal(res$nr.multigraphs,462)
})
