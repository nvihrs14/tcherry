context("weight_junction_tree")
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

tch <- tcherry_step(data, smooth = 0.001)

test_that("weight is correct", {
  expect_equal(weight_junction_tree(tch$cliques, tch$separators, data,
                                    smooth = 0.001),
               tch$weight)
})

cliques <- tch$cliques
cliques[[1]][1] <- "va2"

test_that("error messages work", {
  expect_error(weight_junction_tree(tch$cliques, tch$separators, c(1, 2)),
               "data must be a data frame or a matrix.")
  expect_error(weight_junction_tree(cliques, tch$separators, data),
               "All names in cliques and separators must be column names of data.")
  expect_error(weight_junction_tree(tch$cliques, tch$separators[-1], data),
               "The number of separators must be one less than the number of cliques.
         If a separator is used more than once it should be repeated in the list.")

})
