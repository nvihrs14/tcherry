context("tcherry_complete_search")
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
  expect_error(tcherry_complete_search(data_numeric, 3, smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(tcherry_complete_search(vec, k, smooth = 0.001),
               "data must be a data frame or a matrix.")
  expect_error(tcherry_complete_search(data, vec, smooth = 0.001),
               "k must be a single positive integer.")
  expect_error(tcherry_complete_search(data, 1.1, smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(tcherry_complete_search(data, - 1, smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(tcherry_complete_search(data, 1, smooth = 0.001),
               "k must be a positive integer and at least 2.")
})

tch4 <- tcherry_complete_search(data, 4, smooth = 0.1)
tch5 <- tcherry_complete_search(data, 5, smooth = 0.1)

tch4m <- tcherry_complete_search(data_matrix, 4, smooth = 0.1)
tch5m <- tcherry_complete_search(data_matrix, 5, smooth = 0.1)

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(tcherry_complete_search(data_na, 5, smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("number of models and edges are correct", {
  expect_equal(tch4$n_models, 5915)
  expect_equal(tch5$n_models, 455)
  expect_equal(tch4$model$n_edges, 15)
  expect_equal(tch5$model$n_edges, 18)
  expect_equal(tch4m$n_models, 5915)
  expect_equal(tch5m$n_models, 455)
  expect_equal(tch4m$model$n_edges, 15)
  expect_equal(tch5m$model$n_edges, 18)
})