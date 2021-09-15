test_that("get_degree_seq works", {
  A <-  matrix(c(0, 1, 2,
                 1, 2, 1,
                 2, 1, 2), nrow=3, ncol=3)

  deg <- get_degree_seq(adj = A, type = 'graph')
  expect_equal(deg, c(3,6,7))

  deg <- get_degree_seq(adj = A, type = 'multigraph')
  expect_equal(deg, c(3,4,5))

  expect_error(get_degree_seq(adj = A, type = 'directed'))
})
