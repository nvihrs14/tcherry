# https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17494
suppressWarnings(RNGversion("3.5.0"))

context("cond_independence_test")
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

data_mat <- as.matrix(data)

data_numeric <- data
data_numeric$var2 <- as.numeric(data_numeric$var2)

test_that("error messages work", {
  expect_error(cond_independence_test("var1", "var2",
                                      data = 1:2, smooth = 0.1),
               "data must be a data frame or a matrix.")
  expect_error(cond_independence_test("var1", "var2",
                                      data = data_numeric, smooth = 0.1),
               "Some columns are not characters or factors.")
  expect_error(cond_independence_test("var1", "var2", cond = matrix(1, nrow = 1),
                                      data = data),
               paste("cond must either be a single variable name or a vector,",
                     "possibly empty, of names.", sep = " "))
  expect_error(cond_independence_test("var1", "var2", cond = list("var3"),
                                      data = data),
               paste("cond must either be a single variable name or a vector,",
                     "possibly empty, of names.", sep = " "))  
  expect_error(cond_independence_test(var1 = c("var1", "var3"), var2 = "var2",
                                      data = data),
               "var1 and var2 must each be a single name.")
  expect_error(cond_independence_test("va1", "var2",
                                      data = data, smooth = 0.1),
               "var1, var2 and the variables in cond must be variable names in data.")
  expect_error(cond_independence_test(1, "var2",
                                      data = data, smooth = 0.1),
               "var1, var2 and the variables in cond must be variable names in data.")
  expect_error(cond_independence_test("var1", "var2",
                                      cond = c("v4", "var3"),
                                      data = data, smooth = 0.1),
               "var1, var2 and the variables in cond must be variable names in data.")
  expect_error(cond_independence_test("var1", "var2",
                                      data = data, smooth = 1:2),
               "smooth must be a single non-negative value.")
  expect_error(cond_independence_test("var1", "var2",
                                      data = data, smooth = "a"),
               "smooth must be numeric.")
  expect_error(cond_independence_test("var1", "var2",
                                      data = data, smooth = - 1),
               "smooth must be a non-negative numeric value.")
  expect_error(cond_independence_test("var1", "var2",
                                      data = data),
               "Some probabilities are zero. Consider using the smooth argument.")
})

object <- cond_independence_test("var1", "var4",
                                 data = data, smooth = 0.1)
object_cond <- cond_independence_test("var2", "var3", cond = "var1",
                                      data = data, smooth = 0.1)

object2 <- cond_independence_test("var1", "var4",
                                 data = data_mat, smooth = 0.1)
object_cond2 <- cond_independence_test("var2", "var3", cond = "var1",
                                      data = data_mat, smooth = 0.1)

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(cond_independence_test("var1", "var4", data = data_na,
                                        smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are correct", {
  expect_equal(object$chi_sq_statistic, 0.7512572)
  expect_equal(object$df, 1)
  expect_equal(object_cond$chi_sq_statistic, 0.6297737, tolerance = 1e-7)
  expect_equal(object_cond$df, 8)

  expect_equal(object2$chi_sq_statistic, 0.7512572)
  expect_equal(object2$df, 1)
  expect_equal(object_cond2$chi_sq_statistic, 0.6297737, tolerance = 1e-7)
  expect_equal(object_cond2$df, 8)
})
