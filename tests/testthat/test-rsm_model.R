test_that("rsm_model works", {
  D <- c(2,2,3,3)
  mod1 <- rsm_model(D)
  expect_equal(dim(mod1$m.seq),c(25,10))
  expect_equal(dim(mod1$prob.dists),c(25,4))
  expect_equal(mod1$M$M1[1],0.889)
})
