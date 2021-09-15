test_that("get_edge_multip_seq works", {
  A <-  matrix(c(0, 1, 2,
                 1, 2, 1,
                 2, 1, 2), nrow=3, ncol=3)
  deg <- get_degree_seq(A, 'multigraph')
  expect_equal(dim(get_edge_multip_seq(deg)),c(12,6))
})
