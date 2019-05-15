context("k_tcherry_p_lookahead")
library(tcherry)

set.seed(43)
var1 <- c(sample(c(1, 2), 100, replace = TRUE))
var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
                        prob = c(0.9, 0.1)))
var4 <- c(sample(c(1, 2), 100, replace = TRUE))
var5 <- var2 + var3
var6 <- var1 - var4 + c(sample(c(1, 2), 100, replace = TRUE))

data <- data.frame("var1" = as.character(var1),
                   "var2" = as.character(var2),
                   "var3" = as.character(var3),
                   "var4" = as.character(var4),
                   "var5" = as.character(var5),
                   "var6" = as.character(var6))

data_matrix <- as.matrix(data)

data_numeric <- data
data_numeric[, 3] <- as.numeric(data_numeric[, 3])

vec <- rep(1:2, 5)

test_that("error messages work", {
  expect_error(k_tcherry_p_lookahead(data_numeric, k = 3, p = 3,
                                     smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(k_tcherry_p_lookahead(vec, k = 3, p = 3,
                                     smooth = 0.001),
               "data must be a data frame or a matrix.")
  expect_error(k_tcherry_p_lookahead(data, k = vec, p = 3,
                                     smooth = 0.001),
               "k and p must be single positive integers.")
  expect_error(k_tcherry_p_lookahead(data, k = 1.1, p = 3,
                                     smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(k_tcherry_p_lookahead(data, k= -1, p = 3,
                                     smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(k_tcherry_p_lookahead(data, k = 1, p = 3,
                                     smooth = 0.001),
               "k must be a positive integer and at least 2.")
  expect_error(k_tcherry_p_lookahead(data, k = 3, p = vec,
                                     smooth = 0.001),
               "k and p must be single positive integers.")
  expect_error(k_tcherry_p_lookahead(data, k = 3, p = 1.1,
                                     smooth = 0.001),
               "p must be a positive integer.")
  expect_error(k_tcherry_p_lookahead(data, k= 3, p = - 1,
                                     smooth = 0.001),
               "p must be a positive integer.")
})

tch3_step <- k_tcherry_step(data, 3, smooth = 0.1)
tch3_complete <- tcherry_complete_search(data, 3, smooth = 0.1)$model

tch3_p1 <- k_tcherry_p_lookahead(data, 3, 1, smooth = 0.1)
tch3_pall <- k_tcherry_p_lookahead(data, 3, 4, smooth = 0.1)

tch3_p1m <- k_tcherry_p_lookahead(data_matrix, 3, 1, smooth = 0.1)
tch3_pallm <- k_tcherry_p_lookahead(data_matrix, 3, 4, smooth = 0.1)

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(k_tcherry_p_lookahead(data_na, 3, 1, smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are correct", {
  expect_equal(tch3_p1$adj_matrix, tch3_step$adj_matrix)
  expect_equal(tch3_pall$adj_matrix, tch3_complete$adj_matrix)
  expect_equal(tch3_p1$weight, tch3_step$weight)
  expect_equal(tch3_pall$weight, tch3_complete$weight)
  expect_equal(tch3_p1$n_edges, 9)
  expect_equal(tch3_pall$n_edges, 9)

  expect_equal(tch3_p1m$adj_matrix, tch3_step$adj_matrix)
  expect_equal(tch3_pallm$adj_matrix, tch3_complete$adj_matrix)
  expect_equal(tch3_p1m$weight, tch3_step$weight)
  expect_equal(tch3_pallm$weight, tch3_complete$weight)
  expect_equal(tch3_p1m$n_edges, 9)
  expect_equal(tch3_pallm$n_edges, 9)
})

test_that("Message works", {
  expect_message(k_tcherry_p_lookahead(data, 5, 2, smooth = 0.1),
            paste("Note that p is greather than or equal to the number of cliques",
                  "in the graph. Therfore, all cliques will be added in one step.",
                  sep = " "))
  expect_message(k_tcherry_p_lookahead(data, 5, 3, smooth = 0.1),
            paste("Note that p is greather than or equal to the number of cliques",
                  "in the graph. Therfore, all cliques will be added in one step.",
                  sep = " "))
})