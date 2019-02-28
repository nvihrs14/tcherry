context("tcherry_step")
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

mat_res <- matrix(c(0, 1, 1, 0, 0, 0, 0,
                    1, 0, 1, 0, 1, 1, 0,
                    1, 1, 0, 1, 1, 1, 1,
                    0, 0, 1, 0, 0, 1, 1,
                    0, 1, 1, 0, 0, 0, 0,
                    0, 1, 1, 1, 0, 0, 0,
                    0, 0, 1, 1, 0, 0, 0),
                  nrow = 7)
rownames(mat_res) <- colnames(mat_res) <- names(data)

test_that("results are correct", {
  expect_equal(tcherry_step(data, smooth = 0.001)$adj_matrix,
               mat_res)
  expect_equal(tcherry_step(data, smooth = 0.001)$score, 4.269289,
               tolerance = 1e-6)
})

data_numeric <- data
data_numeric[, 3] <- as.numeric(data_numeric[, 3])

vec <- rep(1:2, 5)

test_that("error messages work", {
  expect_error(tcherry_step(data_numeric, smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(tcherry_step(vec, smooth = 0.001),
               "data must be a data frame or a matrix.")
})
