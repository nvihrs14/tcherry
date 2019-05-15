context("MIk")
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
  expect_error(MIk(c("var1"), vec),
               "data must be a data frame or a matrix.")
  expect_error(MIk(c("va1", "var2"), data),
               "All names in variables must be column names of data.")
  expect_error(MIk(c("var1", "var2", "var3"), data_numeric),
               "Input must be either characters or factors.")
  expect_error(MIk(c("var1", "var2"), data, smooth = c(1, 2)),
               "smooth must be a single non-negative value.")
  expect_error(MIk(c("var1", "var2"), data, smooth = "C"),
               "smooth must be numeric.")
  expect_error(MIk(c("var1", "var2"), data, smooth = - 1),
               "smooth must be a non-negative numeric value.")
  expect_error(MIk(c("var1", "var2", "var3"), data),
               "Some probabilities are zero and therefore MI cannot be calculated.
         Consider using the smooth argument.")
})

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(MIk(c("var1", "var2"), data_na, smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are correct", {
  expect_equal(MIk(c("var1", "var2"), data, smooth = 0.001),
               MI2(data$var1, data$var2, smooth = 0.001))
  expect_equal(MIk(c("var1", "var2", "var3"), data, smooth = 0.001),
               MI3(data$var1, data$var2, data$var3, smooth = 0.001))
  expect_equal(MIk(c("var1", "var2"), data_matrix, smooth = 0.001),
               MI2(data$var1, data$var2, smooth = 0.001))
  expect_equal(MIk(c("var1", "var2", "var3"), data_matrix, smooth = 0.001),
               MI3(data$var1, data$var2, data$var3, smooth = 0.001))
})
