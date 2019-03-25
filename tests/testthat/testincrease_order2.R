context("increase_order2")
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

ChowLiu_cliques <- list(c("var1", "var5"),
                        c("var2", "var5"),
                        c("var3", "var5"),
                        c("var3", "var7"),
                        c("var4", "var6"),
                        c("var5", "var6"))



tcherry_ord_inc <- increase_order2(ChowLiu_cliques, data,
                                              smooth = 0.1)
adj_mat <- matrix(c(0, 0, 1, 1, 1, 1, 0,
                    0, 0, 1, 0, 1, 0, 1,
                    1, 1, 0, 0, 1, 0, 1,
                    1, 0, 0, 0, 0, 1, 0,
                    1, 1, 1, 0, 0, 1, 0,
                    1, 0, 0, 1, 1, 0, 0,
                    0, 1, 1, 0, 0, 0, 0), nrow = 7)
colnames(adj_mat) <- rownames(adj_mat) <- names(data)

cliques <- list(c("var2", "var3", "var5"),
                c("var3", "var5", "var1"),
                c("var5", "var1", "var6"),
                c("var1", "var6", "var4"),
                c("var2", "var3", "var7"))

seps <- list(c("var3", "var5"),
             c("var5", "var1"),
             c("var1", "var6"),
             c("var2", "var3"))


test_that("results are corrects", {
  expect_equal(tcherry_ord_inc$adj_matrix, adj_mat)
  expect_true(compare::compare(tcherry_ord_inc$cliques, cliques,
                               ignoreOrder = TRUE)$result)
  expect_true(compare::compare(tcherry_ord_inc$separators, seps,
                               ignoreOrder = TRUE)$result)
})

data_numeric <- data
data_numeric[, 3] <- as.numeric(data_numeric[, 3])

vec <- rep(1:2, 5)

cliques_small <- list(c("var1", "var2"))

cliques_error <- list(c("var1", "var5"),
                      c("va2", "var5"),
                      c("var3", "var5"),
                      c("var3", "var7"),
                      c("var4", "var6"),
                      c("var5", "var6"))

cliques_error2 <- list(c("var1", "var5"),
                       c("var2", "var5"),
                       c("var3", "var5"),
                       c("var3", "var7"),
                       c("var4"),
                       c("var5", "var6"))

cliques_not_triang <- list(c("var1", "var2"),
                           c("var1", "var3"),
                           c("var2", "var4"),
                           c("var3", "var4"),
                           c("var5", "var7"),
                           c("var6", "var7"))

cliques_wrong_nedges <- list(c("var1", "var3", "var5"),
                             c("var1", "var2", "var4"))

test_that("error messages work", {
  expect_error(increase_order2(ChowLiu_cliques, data_numeric,
                                          smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(increase_order2(ChowLiu_cliques, vec,
                                          smooth = 0.001),
               "data must be a data frame or a matrix.")
  expect_error(increase_order2(vec, data, smooth = 0.1),
               paste("Cliques must be given in a list, each entry containing",
                     "a vector with the names of the variables in the clique.",
                     collapse = " "))
  expect_error(increase_order2(ChowLiu_cliques[-1], data,
                                          smooth = 0.1),
               paste("The column names of data must be the same as the",
                     "variable names in tch_cliq. All variables in data must",
                     "be in at least one clique.", collapse = " "))
  expect_error(increase_order2(cliques_error, data,
                                          smooth = 0.1),
               paste("The column names of data must be the same as the",
                     "variable names in tch_cliq. All variables in data must",
                     "be in at least one clique.", collapse = " "))
  expect_error(increase_order2(cliques_error2, data,
                                          smooth = 0.1),
               paste("tch_cliq should be the cliques of a k'th order t-cherry",
                     "tree. Therefore they should all have the same length k.",
                     collapse = " "))
  expect_error(increase_order2(cliques_small, data[, 1:2],
                                          smooth = 0.1),
               "It takes at least k plus 1 variables to fit a k plus 1'th order t-cherry tree.")
  expect_error(increase_order2(cliques_not_triang, data,
                                          smooth = 0.1),
               paste("The cliques do not come from a triangulated graph.",
                     "The cliques should correspond to a k'th order t-cherry",
                     "tree so it must be triangulated.", collapse = " "))
  expect_error(increase_order2(cliques_wrong_nedges, data[, 1:5],
                                          smooth = 0.1),
               paste("The graph corresponding to the cliques does not have",
                     "the correct number of edges for a k'th order t-cherry",
                     "tree.", collapse = " "))
})
