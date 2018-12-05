context("tcherry")
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

tcherry_tree <- tcherry(data, smooth = 0.1)

tcherry_cliques <- list(c("var2", "var3", "var5"),
                        c("var1", "var2", "var5"),
                        c("var2", "var5", "var6"),
                        c("var3", "var5", "var7"),
                        c("var4", "var5", "var6"))

tcherry_separators <- list(c("var2", "var5"),
                           c("var2", "var5"),
                           c("var3", "var5"),
                           c("var5", "var6"))

test_that("results are corrects", {
  expect_equal(tcherry_tree$adj_tcherry, adj_matrix_tcherry)
  expect_equal(tcherry_tree$cliques, tcherry_cliques)
  expect_equal(tcherry_tree$separators, tcherry_separators)
  })
