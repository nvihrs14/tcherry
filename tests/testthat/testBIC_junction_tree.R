context("BIC_junction_tree")
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
  expect_error(BIC_junction_tree(cliques, separators, data = 1:2),
               "data must be a data frame or a matrix.")
  expect_error(BIC_junction_tree(cliques, separators, data = data_numeric),
               "Some columns are not characters or factors.")
  expect_error(BIC_junction_tree(1:2, separators, data = data),
               paste("Cliques must be given in a list, each entry containing",
                     "a vector with the names of the variables in the clique.",
                     collapse = " "))
  expect_error(BIC_junction_tree(cliques, 1:2, data = data),
               paste("Separators must be given in a list, each entry containing",
                     "a vector with the names of the variables in the separator.",
                     collapse = " "))
  expect_error(BIC_junction_tree(cliques_wrong, separators, data = data),
               paste("The column names of data must be the same as the",
                     "variable names in tch_cliq. All variables in data must",
                     "be in at least one clique.", collapse = " "))
  expect_error(BIC_junction_tree(cliques, separators_wrong, data = data),
               "All variable names in separators should be in data.")
})

test_that("results are correct", {
  expect_equal(BIC_junction_tree(cliques, separators, data, smooth = 0.1),
               -1042.296, tolerance = 10e-5)
})
