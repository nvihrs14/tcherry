# https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17494
suppressWarnings(RNGversion("3.5.0"))

context("increase_order_complete_search")
library(tcherry)

set.seed(43)
var1 <- c(sample(c(1, 2), 100, replace = TRUE))
var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
                        prob = c(0.9, 0.1)))
var4 <- c(sample(c(1, 2), 100, replace = TRUE))

data <- data.frame("var1" = as.character(var1),
                   "var2" = as.character(var2),
                   "var3" = as.character(var3),
                   "var4" = as.character(var4))

data_mat <- as.matrix(data)

ChowLiu_cliques <- list(c("var1", "var2"),
                        c("var1", "var3"),
                        c("var3", "var4"))

tch_complete <- increase_order_complete_search(ChowLiu_cliques, data,
                                              smooth = 0.1)
tch_complete2 <- increase_order_complete_search(ChowLiu_cliques,
                                                data_mat,
                                               smooth = 0.1)

adj_mat <- matrix(c(0, 1, 1, 1,
                    1, 0, 0, 1,
                    1, 0, 0, 1,
                    1, 1, 1, 0), nrow = 4)

colnames(adj_mat) <- rownames(adj_mat) <- names(data)

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(increase_order_complete_search(ChowLiu_cliques, data_na,
                                                smooth = 0.1),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("results are corrects", {
  expect_equal(tch_complete$model$adj_matrix, adj_mat)
  expect_equal(tch_complete$model$n_edges, 5)
  expect_equal(tch_complete$n_models, 3)

  expect_equal(tch_complete2$model$adj_matrix, adj_mat)
  expect_equal(tch_complete2$model$n_edges, 5)
  expect_equal(tch_complete2$n_models, 3)
})

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
  expect_error(increase_order_complete_search(ChowLiu_cliques,
                                              data_numeric,
                                              smooth = 0.001),
               "Some columns are not characters or factors.")
  expect_error(increase_order_complete_search(ChowLiu_cliques, vec,
                                          smooth = 0.001),
               "data must be a data frame or a matrix.")
  expect_error(increase_order_complete_search(vec, data, smooth = 0.1),
               paste("Cliques must be given in a list, each entry",
                     "containing a vector with the names of the",
                     "variables in the clique.",
                     collapse = " "))
  expect_error(increase_order_complete_search(ChowLiu_cliques[- 1], data,
                                              smooth = 0.1),
               paste("The column names of data must be the same as the",
                     "variable names in tch_cliq. All variables in",
                     "data must be in at least one clique.",
                     collapse = " "))
  expect_error(increase_order_complete_search(cliques_error, data,
                                              smooth = 0.1),
               paste("The column names of data must be the same as the",
                     "variable names in tch_cliq. All variables in",
                     "data must be in at least one clique.",
                     collapse = " "))
  expect_error(increase_order_complete_search(cliques_error2, data,
                                              smooth = 0.1),
               paste("tch_cliq should be the cliques of a k'th order",
                     "t-cherry tree. Therefore they should all have",
                     "the same length k.", collapse = " "))
  expect_error(increase_order_complete_search(cliques_small, data[, 1:2],
                                              smooth = 0.1),
               "It takes at least k plus 1 variables to fit a k plus 1'th order t-cherry tree.")
  expect_error(increase_order_complete_search(cliques_not_triang, data,
                                              smooth = 0.1),
               paste("The cliques do not come from a triangulated graph.",
                     "The cliques should correspond to a k'th order t-cherry",
                     "tree so it must be triangulated.", collapse = " "))
  expect_error(increase_order_complete_search(cliques_wrong_nedges, data[, 1:5],
                                              smooth = 0.1),
               paste("The graph corresponding to the cliques does not have",
                     "the correct number of edges for a k'th order t-cherry",
                     "tree.", collapse = " "))
})
