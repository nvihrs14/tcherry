context("diff_edges")
library(tcherry)

m1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
m2_same <- m1
m2_diff <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3)

test_that("number of different edges are correct", {
  expect_equal(diff_edges(m1, m2_same), 0)
  expect_equal(diff_edges(m1, m2_diff), 1)
})
