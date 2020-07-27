# https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17494
suppressWarnings(RNGversion("3.5.0"))

context("tcherry_CL")
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

adj_matrix_tcherry <- matrix(c(0, 1, 0, 0, 1, 0, 0,
                               1, 0, 1, 0, 1, 1, 0,
                               0, 1, 0, 0, 1, 0, 1,
                               0, 0, 0, 0, 1, 1, 0,
                               1, 1, 1, 1, 0, 1, 1,
                               0, 1, 0, 1, 1, 0, 0,
                               0, 0, 1, 0, 1, 0, 0),
                             nrow = 7)

colnames(adj_matrix_tcherry) <- rownames(adj_matrix_tcherry) <-
  names(data)

tcherry_tree <- tcherry_CL(data, smooth = 0.1)
tcherry_tree_matrix <- tcherry_CL(data_matrix, smooth = 0.1)

tcherry_cliques <- list(c("var2", "var3", "var5"),
                        c("var1", "var2", "var5"),
                        c("var2", "var5", "var6"),
                        c("var3", "var5", "var7"),
                        c("var4", "var5", "var6"))

tcherry_separators <- list(c("var2", "var5"),
                           c("var2", "var5"),
                           c("var3", "var5"),
                           c("var5", "var6"))

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(tcherry_CL(data_na, smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are corrects", {
  expect_equal(tcherry_tree$adj_matrix, adj_matrix_tcherry)
  expect_equal(tcherry_tree$cliques, tcherry_cliques)
  expect_equal(tcherry_tree$separators, tcherry_separators)
  expect_equal(tcherry_tree_matrix$adj_matrix, adj_matrix_tcherry)
  expect_equal(tcherry_tree_matrix$cliques, tcherry_cliques)
  expect_equal(tcherry_tree_matrix$separators, tcherry_separators)
  })

data_numeric <- data
data_numeric[, 3] <- as.numeric(data_numeric[, 3])

vec <- rep(1:2, 5)

test_that("error messages work", {
  expect_error(tcherry_CL(data_numeric, smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(tcherry_CL(vec, smooth = 0.001),
               "data must be a data frame or a matrix.")
})
