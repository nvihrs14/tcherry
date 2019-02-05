context("diff_edges")
library(tcherry)

l <- list(1, 2, 3)
letters <- matrix(letters[1:4], 2, 2)
m_notadj <- matrix(c(0, 2, 0, 1, 0, 1, 0, 2, 0), nrow = 3, ncol = 3)
m_notsym <- matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3)
m_loop <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 1), nrow = 3, ncol = 3)

test_that("error messages work", {
  expect_error(diff_edges(l, l), "Arguments must be matrices.")
  expect_error(diff_edges(letters, letters), "Arguments must be numeric.")
  expect_error(diff_edges(m_notadj, m_notadj),
               "Arguments must be adjacency matrices for unweighted graphs.
         Therefore all entries must be 0 or 1.")
  expect_error(diff_edges(m_notsym, m_notsym),
               "Only undirected graphs are supported so arguments must be
         symmetric.")
  expect_error(diff_edges(m_loop, m_loop),
               "Loops are not supported so diagonal must be all 0.")
})

m1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
m2_same <- m1
m2_diff <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3)

test_that("number of different edges are correct", {
  expect_equal(diff_edges(m1, m2_same), 0)
  expect_equal(diff_edges(m1, m2_diff), 1)
})
