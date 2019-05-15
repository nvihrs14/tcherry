context("k_tcherry_step")
library(tcherry)

set.seed(43)
var1 <- c(sample(c(1, 2), 100, replace = TRUE))
var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
                        prob = c(0.9, 0.1)))
var4 <- c(sample(c(1, 2), 100, replace = TRUE))
var5 <- var2 + var3
var6 <- var1 - var4 + c(sample(c(1, 2), 100, replace = TRUE))
var7 <- c(sample(c(1, 2), 100, replace = TRUE))

data <- data.frame("var1" = as.character(var1),
                   "var2" = as.character(var2),
                   "var3" = as.character(var3),
                   "var4" = as.character(var4),
                   "var5" = as.character(var5),
                   "var6" = as.character(var6),
                   "var7" = as.character(var7))

data_matrix <- as.matrix(data)

data_numeric <- data
data_numeric[, 3] <- as.numeric(data_numeric[, 3])

vec <- rep(1:2, 5)

test_that("error messages work", {
  expect_error(k_tcherry_step(data_numeric, 3, smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(k_tcherry_step(vec, k, smooth = 0.001),
               "data must be a data frame or a matrix.")
  expect_error(k_tcherry_step(data, vec, smooth = 0.001),
               "k must be a single positive integer.")
  expect_error(k_tcherry_step(data, 1.1, smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(k_tcherry_step(data, - 1, smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(k_tcherry_step(data, 1, smooth = 0.001),
               "k must be a positive integer and at least 2.")
})

tch2 <- k_tcherry_step(data, 2, smooth = 0.1)
tch3 <- k_tcherry_step(data, 3, smooth = 0.1)

tch2m <- k_tcherry_step(data_matrix, 2, smooth = 0.1)
tch3m <- k_tcherry_step(data_matrix, 3, smooth = 0.1)

CL <- ChowLiu(data, CPTs = FALSE, smooth = 0.1)
tchstep <- tcherry_step(data, smooth = 0.1)

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(k_tcherry_step(data_na, 3, smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are correct", {
  expect_equal(tch2$adj_matrix, CL$skeleton_adj)
  expect_equal(tch3$adj_matrix, tchstep$adj_matrix)
  expect_equal(tch3$n_edges, 11)
  expect_equal(tch2$n_edges, 6)
  expect_equal(tch3$weight, tchstep$weight)
  expect_true(compare::compare(tch3$cliques, tchstep$cliques,
                               ignoreOrder = TRUE)$result)
  expect_true(compare::compare(tch3$separators, tchstep$separators,
                               ignoreOrder = TRUE)$result)

  expect_equal(tch2m$adj_matrix, CL$skeleton_adj)
  expect_equal(tch3m$adj_matrix, tchstep$adj_matrix)
  expect_equal(tch3m$n_edges, 11)
  expect_equal(tch2m$n_edges, 6)
  expect_equal(tch3m$weight, tchstep$weight)
  expect_true(compare::compare(tch3m$cliques, tchstep$cliques,
                               ignoreOrder = TRUE)$result)
  expect_true(compare::compare(tch3m$separators, tchstep$separators,
                               ignoreOrder = TRUE)$result)
})
