context("is_acyclic")
library(tcherry)

adj_matrix_cyclic <- matrix(c(0, 1, 1, 1,
                              1, 0, 0, 1,
                              1, 0, 0, 0,
                              1, 1, 0, 0),
                            nrow = 4)

adj_matrix_acyclic <- matrix(c(0, 0, 1, 1,
                              0, 0, 0, 1,
                              1, 0, 0, 0,
                              1, 1, 0, 0),
                            nrow = 4)
adj_matrix_loop <- matrix(c(1, 0, 1, 1,
                            0, 0, 0, 1,
                            1, 0, 0, 0,
                            1, 1, 0, 0),
                          nrow = 4)
adj_matrix_acyclic_disconnected <- matrix(c(0, 0, 1, 1,
                                            0, 0, 0, 0,
                                            1, 0, 0, 0,
                                            1, 0, 0, 0),
                                          nrow = 4)

test_that("is_acyclic is working", {
  expect_true(is_acyclic(adj_matrix_acyclic))
  expect_true(is_acyclic(adj_matrix_acyclic_disconnected))
  expect_false(is_acyclic(adj_matrix_cyclic))
  expect_error(is_acyclic(adj_matrix_loop),
               "The graph represented by the matrix contains loops.")
})

context("ChowLiu")

set.seed(43)
var1 <- c(sample(c(1, 2), 50, replace = TRUE))
var2 <- var1 + c(sample(c(1, 2), 50, replace = TRUE))
var3 <- var1 + c(sample(c(0, 1), 50, replace = TRUE,
                 prob = c(0.9, 0.1)))
var4 <- c(sample(c(1, 2), 50, replace = TRUE))
data <- data.frame("var1" = as.character(var1),
                   "var2" = as.character(var2),
                   "var3" = as.character(var3),
                   "var4" = as.character(var4))

data_matrix <- as.matrix(data)

data_numeric <- data.frame("var1" = as.character(var1),
                           "var2" = as.character(var2),
                           "var3" = as.character(var3),
                           "var4" = var4)

test_that("Input is specified correctly", {
  expect_error(ChowLiu(data_numeric),
               "Some columns are not characters or factors")
  expect_error(ChowLiu(var1),
               "data must be a data frame or a matrix")
})
