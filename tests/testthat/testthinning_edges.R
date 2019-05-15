context("thinning_edges")
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

cliques <- list(c("var1", "var2", "var3"),
                c("var2", "var3", "var5"),
                c("var5", "var6", "var7"),
                c("var1", "var4"),
                c("var2", "var5", "var6"))

separators <- list(c("var1"),
                   c("var2", "var3"),
                   c("var2", "var5"),
                   c("var5", "var6"))

separators_wrong <- list(c("va1"),
                         c("var2", "var3"),
                         c("var2", "var5"),
                         c("var5", "var6"))

cliques_wrong <- list(c("var1", "var2", "var3"),
                      c("var2", "var3", "var5"),
                      c("var5", "var6", "var7"),
                      c("var1", "var4"),
                      c("var2", "var5", "var8"))

test_that("error messages work", {
  expect_error(thinning_edges(cliques, separators, data = 1:2),
               "data must be a data frame or a matrix.")
  expect_error(thinning_edges(cliques, separators, data = data_numeric),
               "Some columns are not characters or factors.")
  expect_error(thinning_edges(1:2, separators, data = data),
               paste("Cliques must be given in a list, each entry containing",
                     "a vector with the names of the variables in the clique.",
                     collapse = " "))
  expect_error(thinning_edges(cliques, 1:2, data = data),
               paste("Separators must be given in a list, each entry containing",
                     "a vector with the names of the variables in the separator.",
                     collapse = " "))
  expect_error(thinning_edges(cliques_wrong, separators, data = data),
               paste("The column names of data must be the same as the",
                     "variable names in tch_cliq. All variables in data must",
                     "be in at least one clique.", collapse = " "))
  expect_error(thinning_edges(cliques, separators_wrong, data = data),
               "All variable names in separators should be in data.")
  expect_error(thinning_edges(cliques, separators, data = data, alpha = 1:2),
               "alpha must be a single non-negative value.")
  expect_error(thinning_edges(cliques, separators, data = data, alpha = "a"),
               "alpha must be numeric.")
  expect_error(thinning_edges(cliques, separators, data = data, alpha = - 1),
               "alpha must be a positive numeric value.")
})

graph <- thinning_edges(cliques, separators, data = data, smooth = 0.1)
graph2 <- thinning_edges(cliques, separators, data = data_mat,
                         smooth = 0.1)

target_mat <- matrix(c(0, 1, 1, 0, 0, 0, 0,
                       1, 0, 1, 0, 1, 0, 0,
                       1, 1, 0, 0, 1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0,
                       0, 1, 1, 0, 0, 1, 0,
                       0, 0, 0, 0, 1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0), nrow = 7)

colnames(target_mat) <- rownames(target_mat) <- names(data)

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(thinning_edges(cliques, separators, data_na, smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are correct", {
  expect_equal(graph$adj_matrix, target_mat)
  expect_equal(graph$n_edges_removed, 4)
  expect_equal(graph$n_edges, 6)
  expect_equal(graph2$adj_matrix, target_mat)
  expect_equal(graph2$n_edges_removed, 4)
  expect_equal(graph2$n_edges, 6)
})