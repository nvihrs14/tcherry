context("diff_edges_tch")
library(tcherry)

l <- list(1, 2, 3)
letters <- matrix(letters[1:4], 2, 2)
m_notadj <- matrix(c(0, 2, 0, 1, 0, 1, 0, 2, 0), nrow = 3, ncol = 3)
m_notsym <- matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3)
m_loop <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 1), nrow = 3, ncol = 3)
m2 <- matrix(c(0, 1, 1, 0), 2, 2)
m3 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)

m_nocol <- m3
rownames(m_nocol) <- letters[1:3]

m_diffnames <- m_nocol
colnames(m_diffnames) <- letters[4:6]

m_order_names <- m_diffnames
colnames(m_order_names) <- c("a", "c", "b")

m_names <- m_nocol
colnames(m_names) <- rownames(m_names)

m_names_2 <- m3
rownames(m_names_2) <- colnames(m_names_2) <- letters[4:6]

m_nonames <- matrix(c(0, 1, 1, 0), 2, 2)

test_that("error messages work", {
  expect_error(diff_edges_tch(l, l), "Arguments must be matrices.")
  expect_error(diff_edges_tch(letters, letters), "Arguments must be numeric.")
  expect_error(diff_edges_tch(m2, m3),
               "The matrices must have the same dimensions.")
  expect_error(diff_edges_tch(m_notadj, m_notadj),
               "Arguments must be adjacency matrices for unweighted graphs.
         Therefore all entries must be 0 or 1.")
  expect_error(diff_edges_tch(m_notsym, m_notsym),
               "Only undirected graphs are supported so arguments must be
         symmetric. This includes that rownames must equal colnames.")
  expect_error(diff_edges_tch(m_loop, m_loop),
               "Loops are not supported so diagonal must be all 0.")

  expect_error(diff_edges_tch(m_nocol, m_nocol),
               "Only undirected graphs are supported so arguments must be
         symmetric. This includes that rownames must equal colnames.")
  expect_error(diff_edges_tch(m_diffnames, m_diffnames),
               "Only undirected graphs are supported so arguments must be
         symmetric. This includes that rownames must equal colnames.")
  expect_error(diff_edges_tch(m_order_names, m_order_names),
               "Only undirected graphs are supported so arguments must be
         symmetric. This includes that rownames must equal colnames.")
  expect_error(diff_edges_tch(m_names, m_names_2),
               "The node names must be the same in both graphs.")
  expect_error(diff_edges_tch(m_nonames, m_nonames),
               "The matrices must be named.")
})

m1 <- matrix(c(0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0),
             nrow = 4, ncol = 4)
m2 <- matrix(c(0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0),
             nrow = 4, ncol = 4)

rownames(m1) <- colnames(m1) <- letters[1:4]
rownames(m2) <- colnames(m2) <- letters[1:4]


m1_names <- m1
colnames(m1_names) <- rownames(m1_names) <- letters[1:4]
m2_names <- m1_names[c("a", "c", "b", "d"), c("a", "c", "b", "d")]

test_that("number of different edges are correct", {
  expect_equal(diff_edges_tch(m1, m1), 0)
  expect_equal(diff_edges_tch(m1, m2), 1)
  expect_equal(diff_edges_tch(m1_names, m2_names), 0)
})
